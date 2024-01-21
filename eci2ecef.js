// Ecef to Eci coordinate frame conversion code
// Javascript codes are implemented from Matlab source codes from
// https://github.com/Spacecraft-Code/Vallado which was imported from celestrak at 2017-06-30
// Note: there are many methods for converting Ecef to Eci in the original source
//       however, this project only implements:
//       * class equinox based, 2000b method for precession-nutation effects (e.g., iau00f2i 'b' option)
// tested with ex3_1415.m in matlab source codes
// ex3_1415.m results are verified with
//  https://hpiers.obspm.fr/eop-pc/index.php?index=rotation&lang=en
//  Matlab - HPIERS errors are up to few milimeters
//  Matlab - JS errors are around 1e-12
//
//
// TODO: import from https://celestrak.org/software/vallado/matlab.zip as of 2023 May 09

const vlib = require("../vectorLib/vector");
const fs = require("fs");
let folderPath = "./";

// test code:

// // download latest eop:

//const t = new Date("2019-01-01T00:00:00.000Z") / 1000.0;
// checkAndUpdateEopFile(t);
// const eop = readEopParameters(t);

// // unix time in seconds:
// // const t = new Date("2014-01-01T00:00:00.000Z") / 1000.0;
// const recef = vlib.vec(-1033.479383, 7901.2952754, 6380.3565958);
// const vecef = vlib.vec(-3.22563652, -2.87245145, 5.531924446);
// const aecef = vlib.vec(0.001, 0.002, 0.003);
// console.log(ecef2eci(t, recef, vecef, aecef));
// // end of test code

// ecef to eci given unix time in seconds (UTC0)
function ecef2eci(tUnixSec, recef, vecef, aecef) {
  // define undefined inputs:
  if (recef == undefined) recef = vlib.vec(0, 0, 0);
  if (vecef == undefined) vecef = vlib.vec(0, 0, 0);
  if (aecef == undefined) aecef = vlib.vec(0, 0, 0);

  const conv = Math.PI / (180 * 3600);
  const dat = eop.dat; // TAI - UTC
  const dut1 = eop.dut1; // UT1 - UTC
  const xp = eop.xp * conv; //  % " to rad
  const yp = eop.yp * conv;
  const lod = eop.lod;
  const jd = unix2jd(tUnixSec);

  // perform time conversions
  const jdut1 = jd + dut1 / 86400.0;
  // const jdtai = jd + dat/86400.0; // jd of atomic time
  const jdtt = jd + (dat + 32.184) / 86400.0; // jd of terrestrial time
  const ttt = (jdtt - 2451545.0) / 36525.0;

  return iau00f2i(recef, vecef, aecef, ttt, jdut1, lod, xp, yp);
}

// %                           function iau00f2i
// %
// %  this function transforms a vector from the earth fixed (itrf) frame, to
// %    the eci mean equator mean equinox (gcrf).
// %
// %  inputs          description                    range / units
// %    recef       - position vector earth fixed    km
// %    vecef       - velocity vector earth fixed    km/s
// %    aecef       - acceleration vector earth fixedkm/s2
// %    ttt         - julian centuries of tt         centuries
// %    jdut1       - julian date of ut1             days from 4713 bc
// %    lod         - excess length of day           sec -> required for v and a conversion
// %    xp          - polar motion coefficient       arc sec
// %    yp          - polar motion coefficient       arc sec
// %
// %  outputs       :
// %    reci        - position vector eci            km
// %    veci        - velocity vector eci            km/s
// %    aeci        - acceleration vector eci        km/s2
// %
// %  locals        :
// %    pm          - transformation matrix for itrf-pef
// %    st          - transformation matrix for pef-ire
// %    nut         - transformation matrix for ire-gcrf
// %
// %  coupling      :
// %   iau00pm      - rotation for polar motion      itrf-pef
// %   iau00era     - rotation for earth rotyation   pef-ire
// %   iau00xys     - rotation for prec/nut          ire-gcrf
//
// %
// %  references    :
// %    vallado       2004, 205-219
function iau00f2i(recef, vecef, aecef, ttt, jdut1, lod, xp, yp) {
  // % ---- class equinox based, 2000b
  const pn = iau00pnb(ttt);

  const st = iau00gst(
    jdut1,
    ttt,
    pn.deltapsi,
    pn.l,
    pn.l1,
    pn.f,
    pn.d,
    pn.omega,
    pn.lonmer,
    pn.lonven,
    pn.lonear,
    pn.lonmar,
    pn.lonjup,
    pn.lonsat,
    pn.lonurn,
    pn.lonnep,
    pn.precrate
  );

  const pm = polarm(xp, yp, ttt);

  //% ---- setup parameters for velocity transformations
  const thetasa = 7.29211514670698e-5 * (1.0 - lod / 86400.0);
  const omegaearth = vlib.vec(0, 0, thetasa);

  // reci = Apef2eci * Aecef2pef*recef
  const Aecef2pef = pm;
  const Apef2eci = vlib.MxM(pn.pnb, st);
  const Aecef2eci = vlib.MxM(Apef2eci, Aecef2pef);
  const rpef = vlib.Mxv(Aecef2pef, recef);
  const reci = vlib.Mxv(Apef2eci, rpef);

  // veci = Apef2eci * (Aecef2pef*vecef + omegaearth x rpef)
  const temp = vlib.cross(omegaearth, rpef);
  const vpef = vlib.Mxv(Aecef2pef, vecef);
  const veci = vlib.Mxv(Apef2eci, vlib.add(vpef, temp));

  // aeci = Apef2eci * (Aecef2pef*aecef + omegaearth x omegaearth x rpef + 2*omegaearth x vpef)
  const apef = vlib.add(
    vlib.Mxv(Aecef2pef, aecef),
    vlib.add(
      vlib.cross(omegaearth, temp),
      vlib.mult(2, vlib.cross(omegaearth, vpef))
    )
  );
  const aeci = vlib.Mxv(Apef2eci, apef);

  return {
    reci: reci,
    veci: veci,
    aeci: aeci,
    Aecef2eci: Aecef2eci,
    Aecef2pef: Aecef2pef,
    omegaearth: omegaearth,
  };
}

// %                           function iau00pnb
// %
// %  this function calulates the transformation matrix that accounts for the
// %    effects of precession-nutation in the iau2000b theory.
// %
// %  inputs          description                    range / units
// %    ttt         - julian centuries of tt
// %
// %  outputs       :
// %    nut         - transformation matrix for ire-gcrf
// %    deltapsi    - change in longitude            rad
// %    l           - delaunay element               rad
// %    ll          - delaunay element               rad
// %    f           - delaunay element               rad
// %    d           - delaunay element               rad
// %    omega       - delaunay element               rad
//
// %  references    :
// %    vallado       2004, 212-214
function iau00pnb(ttt) {
  // % " to rad
  const convrt = Math.PI / (180.0 * 3600.0);
  // const deg2rad = Math.PI / 180.0;

  // const ttt2 = ttt * ttt;
  // const ttt3 = ttt2 * ttt;
  // const ttt4 = ttt2 * ttt2;
  // const ttt5 = ttt3 * ttt2;

  // % obtain data for calculations form the 2000b theory
  const args = fundarg(ttt);
  const pnConst = generateIau00in_constants();
  let pnsum = 0.0;
  let ensum = 0.0;
  let i;
  for (i = 76; i >= 0; i--) {
    const tempval =
      pnConst.apni[i][0] * args.l +
      pnConst.apni[i][1] * args.l1 +
      pnConst.apni[i][2] * args.f +
      pnConst.apni[i][3] * args.d +
      pnConst.apni[i][4] * args.omega;
    pnsum +=
      (pnConst.apn[i][0] + pnConst.apn[i][1] * ttt) * Math.sin(tempval) +
      (pnConst.apn[i][4] + pnConst.apn[i][5] * ttt) * Math.cos(tempval);
    ensum +=
      (pnConst.apn[i][2] + pnConst.apn[i][3] * ttt) * Math.cos(tempval) +
      (pnConst.apn[i][6] + pnConst.apn[i][7] * ttt) * Math.sin(tempval);
  }

  // % ------ form the planetary arguments
  let pplnsum = -0.000135 * convrt; //% " to rad
  let eplnsum = 0.000388 * convrt;

  // %  add planetary and luni-solar components.
  let deltapsi = pnsum + pplnsum;
  let deltaeps = ensum + eplnsum;
  const p = precess(ttt);

  const oblo = 84381.406 * convrt; //% " to rad

  // % ----------------- find nutation matrix ----------------------
  // % mean to true
  const a1 = vlib.rot1mat(p.ea + deltaeps);
  const a2 = vlib.rot3mat(deltapsi);
  const a3 = vlib.rot1mat(-p.ea);

  // % j2000 to date (precession)
  const a4 = vlib.rot3mat(-p.xa);
  const a5 = vlib.rot1mat(p.wa);
  const a6 = vlib.rot3mat(p.psia);
  const a7 = vlib.rot1mat(-oblo);

  // % icrs to j2000
  const a8 = vlib.rot1mat(-0.0068192 * convrt);
  const a9 = vlib.rot2mat(0.041775 * Math.sin(oblo) * convrt);
  // ?%      a9  = rot2mat(0.0166170*convrt);
  const a10 = vlib.rot3mat(0.0146 * convrt);

  // a10*a9*a8*a7*a6*a5*a4
  const prec = vlib.MxM(
    a10,
    vlib.MxM(a9, vlib.MxM(a8, vlib.MxM(a7, vlib.MxM(a6, vlib.MxM(a5, a4)))))
  );

  // a3*a2*a1
  const nut = vlib.MxM(a3, vlib.MxM(a2, a1));

  const pnb = vlib.MxM(prec, nut);

  return {
    deltapsi: deltapsi,
    pnb: pnb,
    //prec: prec,
    nut: nut,
    l: args.l,
    l1: args.l1,
    f: args.f,
    d: args.d,
    omega: args.omega,
    lonmer: args.lonmer,
    lonven: args.lonven,
    lonear: args.lonear,
    lonmar: args.lonmar,
    lonjup: args.lonjup,
    lonsat: args.lonsat,
    lonurn: args.lonurn,
    lonnep: args.lonnep,
    precrate: args.precrate,
  };
}

// %                           function fundarg
// %
// %  this function calulates the delauany variables and planetary values for
// %  several theories.
// %
// %  author        : david vallado                  719-573-2600   16 jul 2004
// %
// %  revisions
// %    vallado     - consolidate with iau 2000                     14 feb 2005
// %
// %  inputs          description                    range / units
// %    ttt         - julian centuries of tt
// %    opt         - method option                  '10', '02', '96', '80'
// %
// %  outputs       :
// %    l           - delaunay element               rad
// %    ll          - delaunay element               rad
// %    f           - delaunay element               rad
// %    d           - delaunay element               rad
// %    omega       - delaunay element               rad
// %    planetary longitudes                         rad
// %
// %  locals        :
// %    ttt2,ttt3,  - powers of ttt
// %
// %  coupling      :
// %    none        -
// %
// %  references    :
// %    vallado       2004, 212-214
// %
function fundarg(ttt) {
  let deg2rad = Math.PI / 180.0;

  // iau 2000b theory:
  // % ------ form the delaunay fundamental arguments in deg
  let l = 134.96340251 + (1717915923.2178 * ttt) / 3600.0;
  let l1 = 357.52910918 + (129596581.0481 * ttt) / 3600.0;
  let f = 93.27209062 + (1739527262.8478 * ttt) / 3600.0;
  let d = 297.85019547 + (1602961601.209 * ttt) / 3600.0;
  let omega = 125.04455501 + (-6962890.5431 * ttt) / 3600.0;

  // % ------ form the planetary arguments in deg
  let lonmer = 0.0;
  let lonven = 0.0;
  let lonear = 0.0;
  let lonmar = 0.0;
  let lonjup = 0.0;
  let lonsat = 0.0;
  let lonurn = 0.0;
  let lonnep = 0.0;
  let precrate = 0.0;

  // % ---- convert units to rad
  l = (l % 360.0) * deg2rad; // % rad
  l1 = (l1 % 360.0) * deg2rad;
  f = (f % 360.0) * deg2rad;
  d = (d % 360.0) * deg2rad;
  omega = (omega % 360.0) * deg2rad;

  lonmer = (lonmer % 360.0) * deg2rad; //  % rad
  lonven = (lonven % 360.0) * deg2rad;
  lonear = (lonear % 360.0) * deg2rad;
  lonmar = (lonmar % 360.0) * deg2rad;
  lonjup = (lonjup % 360.0) * deg2rad;
  lonsat = (lonsat % 360.0) * deg2rad;
  lonurn = (lonurn % 360.0) * deg2rad;
  lonnep = (lonnep % 360.0) * deg2rad;
  precrate = (precrate % 360.0) * deg2rad;

  return {
    // % ------ form the delaunay fundamental arguments in deg
    l: l,
    l1: l1,
    f: f,
    d: d,
    omega: omega,
    // % ------ form the planetary arguments in deg
    lonmer: lonmer,
    lonven: lonven,
    lonear: lonear,
    lonmar: lonmar,
    lonjup: lonjup,
    lonsat: lonsat,
    lonurn: lonurn,
    lonnep: lonnep,
    precrate: precrate,
  };
}

// %                           function iau00in
// %  this function initializes the matricies needed for iau 2000 reduction
// %    calculations. the routine uses the files listed as inputs, but they are
// %    are not input to the routine as they are static files.
// %
// %
// %  this function initializes the matricies needed for iau 2000 reduction
// %    calculations. the routine uses the files listed as inputs, but they are
// %    are not input to the routine as they are static files.
// %  inputs          description                    range / units
// %    iau00s.dat  - file for s coefficient
// %    iau00n.dat  - file for nutation coefficients
// %
// %  outputs       :
// %    ass0        - real coefficients for s        rad
// %    a0si        - integer coefficients for s
// %    apn         - real coefficients for nutation rad
// %    apni        - integer coefficients for nutation
// Generate json for constants:
function generateIau00in_constants() {
  // % " to rad
  const convrtu = (0.000001 * Math.PI) / (180.0 * 3600.0); // % if micro arcsecond
  const convrtm = (0.001 * Math.PI) / (180.0 * 3600.0); // % if milli arcsecond

  const as = file2Matrix(folderPath+"iau00s.dat");
  // multiply 2nd and 3rd columns with convrtu and assign to ass0
  const ass0 = as
    .map((row) => row.slice(1, 3))
    .map((row) => row.map((col) => col * convrtu));
  const a0si = as.map((row) => row.slice(3));

  const an = file2Matrix(folderPath+"iau03n.dat");
  // first 5 columns are assigned to apni
  const apni = an.map((row) => row.slice(0, 5));
  // rest of collumns are assigned to apn
  const apn = an
    .map((row) => row.slice(6))
    .map((row) => row.map((col) => col * convrtm));

  const out = {
    ass0: ass0,
    a0si: a0si,
    apn: apn,
    apni: apni,
  };
  return out;
  // fs.writeFileSync("pnConstants.txt", JSON.stringify(out), "utf8");
}

// %                           function precess
// %
// %  this function calulates the transformation matrix that accounts for the effects
// %    of precession. both the 1980 and 2000 theories are handled. note that the
// %    required parameters differ a little.
// %
// %  inputs          description                    range / units
// %    ttt         - julian centuries of tt
// %    opt         - method option                  '01', '02', '96', '80'
// %
// %  outputs       :
// %    prec        - transformation matrix for mod - j2000 (80 only)
// %    psia        - cannonical precession angle    rad    (00 only)
// %    wa          - cannonical precession angle    rad    (00 only)
// %    ea          - cannonical precession angle    rad    (00 only)
// %    xa          - cannonical precession angle    rad    (00 only)
// %
// %  locals        :
// %    ttt2        - ttt squared
// %    ttt3        - ttt cubed
// %    zeta        - precession angle               rad
// %    z           - precession angle               rad
// %    theta       - precession angle               rad
// %    oblo        - obliquity value at j2000 epoch "%
// %
// %  references    :
// %    vallado       2004, 214-216, 219-221
// %
// % [prec,psia,wa,ea,xa] = precess ( ttt, opt );
function precess(ttt) {
  const convrt = Math.PI / (180.0 * 3600.0);
  // const ttt2 = ttt * ttt;
  // const ttt3 = ttt2 * ttt;

  // % ------------------ iau 03 precession angles -------------------
  const oblo_ = 84381.406;
  const psia =
    convrt *
    (((((-0.0000000951 * ttt + 0.000132851) * ttt - 0.00114045) * ttt -
      1.0790069) *
      ttt +
      5038.481507) *
      ttt);
  const wa =
    convrt *
    (((((0.0000003337 * ttt - 0.000000467) * ttt - 0.00772503) * ttt +
      0.0512623) *
      ttt -
      0.025754) *
      ttt +
      oblo_);
  const ea =
    convrt *
    (((((-0.0000000434 * ttt - 0.000000576) * ttt + 0.0020034) * ttt -
      0.0001831) *
      ttt -
      46.836769) *
      ttt +
      oblo_);
  const xa =
    convrt *
    ((((-0.000000056 * ttt + 0.000170663) * ttt - 0.00121197) * ttt -
      2.3814292) *
      ttt +
      10.556403) *
    ttt;
  // const zeta =  convrt * ((((( - 0.0000003173 * ttt - 0.000005971 ) * ttt + 0.01801828 ) * ttt + 0.2988499 ) * ttt + 2306.083227 ) * ttt + 2.650545);
  // const theta=  convrt * (((( - 0.0000001274 * ttt - 0.000007089 ) * ttt - 0.04182264 ) * ttt - 0.4294934 ) * ttt + 2004.191903 ) * ttt;
  // const z    =  convrt * (((((   0.0000002904 * ttt - 0.000028596 ) * ttt + 0.01826837 ) * ttt + 1.0927348 ) * ttt + 2306.077181 ) * ttt - 2.650545);
  // const oblo =  convrt * oblo_;
  // const a4   =  vlib.rot3mat(-xa);
  // const a5   =  vlib.rot1mat(wa);
  // const a6   =  vlib.rot3mat(psia);
  // const a7   =  vlib.rot1mat(-oblo);
  // const prec =  vlib.MxM(a7, vlib.MxM(a6, vlib.MxM(a5, a4)));
  return {
    // prec: prec,
    psia: psia,
    wa: wa,
    ea: ea,
    xa: xa,
  };
}

// %                           function iau00gst
// %
// %  this function finds the iau2000 greenwich sidereal time.
// %
// %  inputs          description                    range / units
// %    jdut1       - julian date of ut1             days from 4713 bc
// %    ttt         - julian centuries of tt
// %    deltapsi    - change in longitude            rad
// %    l           - delaunay element               rad
// %    ll          - delaunay element               rad
// %    f           - delaunay element               rad
// %    d           - delaunay element               rad
// %    omega       - delaunay element               rad
// %    many others for planetary values             rad
// %
// %  outputs       :
// %    gst         - greenwich sidereal time        0 to twopi rad
// %    st          - transformation matrix
// %
// %  locals        :
// %    temp        - temporary variable for reals   rad
// %    tut1d       - days from the jan 1, 2000 12 h epoch (ut1)
// %
// %  coupling      :
// %    iau00in     - initialize the data arrays
// %
// %  references    :
// %    vallado       2004, 216
// %
// % [gst,st] = iau00gst(jdut1, ttt, deltapsi, l, l1, f, d, omega, ...
// %            lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate);
// % -----------------------------------------------------------------------------
function iau00gst(
  jdut1,
  ttt,
  deltapsi,
  l,
  l1,
  f,
  d,
  omega,
  lonmer,
  lonven,
  lonear,
  lonmar,
  lonjup,
  lonsat,
  lonurn,
  lonnep,
  precrate
) {
  const deg2rad = Math.PI / 180.0;
  // % " to rad
  const convrt = Math.PI / (180.0 * 3600.0);

  // read agst and agsti from iau00gs.dat file:
  const convrtu = (0.000001 * Math.PI) / (180.0 * 3600.0); // % if micro arcsecond
  const ag = file2Matrix(folderPath+"iau00gs.dat");
  const agst = ag
    .map((row) => row.slice(1, 3))
    .map((row) => row.map((col) => col * convrtu));
  const agsti = ag.map((row) => row.slice(3));
  // fs.writeFileSync("pnConstants2.txt", JSON.stringify({agsti:agsti, agst:agst}), "utf8");

  const ttt2 = ttt * ttt;
  const ttt3 = ttt2 * ttt;
  const ttt4 = ttt2 * ttt2;

  // % mean obliquity of the ecliptic
  let epsa = 84381.448 - 46.84024 * ttt - 0.00059 * ttt2 + 0.001813 * ttt3; // % "
  epsa = (epsa / 3600.0) % 360.0; // % deg
  epsa = epsa * deg2rad; // % rad

  // %  evaluate the ee complementary terms
  let gstsum0 = 0.0;
  for (let i = 32; i >= 0; i--) {
    const tempval =
      agsti[i][0] * l +
      agsti[i][1] * l1 +
      agsti[i][2] * f +
      agsti[i][3] * d +
      agsti[i][4] * omega +
      agsti[i][5] * lonmer +
      agsti[i][6] * lonven +
      agsti[i][7] * lonear +
      agsti[i][8] * lonmar +
      agsti[i][9] * lonjup +
      agsti[i][10] * lonsat +
      agsti[i][11] * lonurn +
      agsti[i][12] * lonnep +
      agsti[i][13] * precrate;
    gstsum0 += agst[i][0] * Math.sin(tempval) + agst[i][1] * Math.cos(tempval); // % rad
  }

  // gstsum1 = 0.0;
  // for j = 1: -1 : 1
  //     i = 33 + j;
  let i = 33;
  const tempval =
    agsti[i][0] * l +
    agsti[i][1] * l1 +
    agsti[i][2] * f +
    agsti[i][3] * d +
    agsti[i][4] * omega +
    agsti[i][5] * lonmer +
    agsti[i][6] * lonven +
    agsti[i][7] * lonear +
    agsti[i][8] * lonmar +
    agsti[i][9] * lonjup +
    agsti[i][10] * lonsat +
    agsti[i][11] * lonurn +
    agsti[i][12] * lonnep +
    agsti[i][13] * precrate;
  const gstsum1 =
    (agst[i][0] * Math.sin(tempval) + agst[i][1] * Math.cos(tempval)) * ttt;

  const eect2000 = gstsum0 + gstsum1 * ttt; // % rad

  // % equation of the equinoxes
  const ee2000 = deltapsi * Math.cos(epsa) + eect2000; //  % rad

  // %  earth rotation angle
  const tut1d = jdut1 - 2451545.0;
  const twopi = 2 * Math.PI;
  let era = twopi * (0.779057273264 + 1.00273781191135448 * tut1d);
  era = era % twopi; //  % rad

  // %  greenwich mean sidereal time, iau 2000.
  const gmst2000 =
    era +
    (0.014506 +
      4612.15739966 * ttt +
      1.39667721 * ttt2 -
      0.00009344 * ttt3 +
      0.00001882 * ttt4) *
      convrt; // % " to rad

  const gst = gmst2000 + ee2000; // % rad

  // % transformation matrix
  const st = vlib.rot3mat(-gst);
  return st;
}

// %                           function polarm
// %
// %  this function calulates the transformation matrix that accounts for polar
// %    motion. both the 1980 and 2000 theories are handled. note that the rotation
// %    order is different between 1980 and 2000 .
// %
// %  author        : david vallado                  719-573-2600   25 jun 2002
// %
// %  revisions
// %    vallado     - consolidate with iau 2000                     14 feb 2005
// %
// %  inputs          description                    range / units
// %    xp          - polar motion coefficient       rad
// %    yp          - polar motion coefficient       rad
// %    ttt         - julian centuries of tt (00 theory only)
// %
// %  outputs       :
// %    pm          - transformation matrix for ecef - pef
// %
// %  locals        :
// %    convrt      - conversion from arcsec to rad
// %    sp          - s prime value
// %
// %  references    :
// %    vallado       2004, 207-209, 211, 223-224
// %
// % [pm] = polarm ( xp, yp, ttt, opt );
// % ----------------------------------------------------------------------------
function polarm(xp, yp, ttt) {
  const cosxp = Math.cos(xp);
  const sinxp = Math.sin(xp);
  const cosyp = Math.cos(yp);
  const sinyp = Math.sin(yp);
  const convrt = Math.PI / (3600.0 * 180.0);
  const sp = -47.0e-6 * ttt * convrt;
  const cossp = Math.cos(sp);
  const sinsp = Math.sin(sp);

  // % form the matrix
  return [
    [
      cosxp * cossp,
      -cosyp * sinsp + sinyp * sinxp * cossp,
      -sinyp * sinsp - cosyp * sinxp * cossp,
    ],
    [
      cosxp * sinsp,
      cosyp * cossp + sinyp * sinxp * sinsp,
      sinyp * cossp - cosyp * sinxp * sinsp,
    ],
    [sinxp, -sinyp * cosxp, cosyp * cosxp],
  ];
}

// // supplementary function to parse dat files
function file2Matrix(filename) {
  const lines = readFile(filename);
  const rows = lines.length;
  let out = [];
  for (let row = 0; row < rows; row++)
    // Use a regular expression to split values by any number of spaces,
    // then, remove empty values from the array
    out[row] = lines[row]
      .split(/\s+/)
      .filter((value) => value !== "" && value != "");
  return out
    .map((row) => row.map((col) => parseFloat(col)))
    .filter((row) =>
      row.some(
        (element) => element !== undefined && element !== null && element !== ""
      )
    );
}

// supplementary function to read a file convert to line array
function readFile(filename) {
  const data = fs.readFileSync(filename, {
    encoding: "utf8",
    flag: "r",
  });
  // return lines as array
  return data.split("\n");
}

// unix seconds to julian days
function unix2jd(unix) {
  return unix / 86400.0 + 2440587.5;
}

// unix seconds to modified julian days
function unix2mjd(unix) {
  return unix / 86400.0 + 40587.0;
}

// modified julian days to unix seconds
function mjd2unix(mjd) {
  return (mjd - 40587.0) * 86400.0;
}

// download Eop file from celestrack
async function updateEopFile() {
  const url = "https://celestrak.org/SpaceData/EOP-Last5Years.csv";
  await fetch(url)
    .then(async (response) => {
      let data = await response.text();
      fs.writeFileSync(folderPath+"eopData.csv", data, "utf8");
    })
    .catch((err) =>
      console.log(
        "could not retrieve EOP files from celestrak.org,\n Error Code:\n",
        err
      )
    );
}

// check eop file and update if it is outdated
// tCheck is in unix seconds
function checkAndUpdateEopFile(tCheck) {
  const lines = readFile(folderPath+"eopData.csv");
  const lastLineValues = lines[lines.length - 2].split(",");
  const lastDate = new Date(1000 * mjd2unix(lastLineValues[1])) / 1000.0;
  if (tCheck == undefined) tCheck = new Date() / 1000.0;
  const isEopOutdated = lastDate < tCheck;
  if (isEopOutdated) {
    console.log("EOP file is upating");
    updateEopFile();
  }
}

// read EOP parameters given unix seconds
// parameter at begin of the day of the unix seconds is used,
// no data interpolation is made
function readEopParameters(tUnix) {
  const lines = readFile(folderPath+"eopData.csv");
  const firstLineValues = lines[1].split(",");
  const lastLineValues = lines[lines.length - 3].split(",");
  const mjdFirst = new Date(parseFloat(firstLineValues[1]));
  const mjdLast = new Date(parseFloat(lastLineValues[1]));
  const mjd = unix2mjd(tUnix);
  if (mjd < mjdFirst || mjd > mjdLast) {
    console.log(
      "Error: requested time is earlier than available EOP file\n Earliest Time in EOP: ",
      new Date(1000 * mjd2unix(mjdFirst)),
      "\n requested time: ",
      new Date(1000 * mjd2unix(mjd)),
      "\n Using default parameters"
    );
    return {
      xp: 0.0,
      yp: 0.0,
      dut1: 0.0, // UT1 - UTC
      lod: 0.0, // length of day difference
      dat: 37, // TAI - UTC
    };
  } else {
    const diff = Math.floor(mjd) - mjdFirst;
    const v = lines[1 + diff].split(",");
    return {
      xp: parseFloat(v[2]),
      yp: parseFloat(v[3]),
      dut1: parseFloat(v[4]), // UT1 - UTC
      lod: parseFloat(v[5]), // length of day difference
      dat: parseFloat(v[10]), // TAI - UTC
    };
  }
}

/** Class for applying ECEF to ECI or vice-versa coordinate transforms.
 * Retrieving Earth Orientation parameters from web
 * @class
 */
class EarthOrientation {
  // public members:
  /**
   * @public
   * @member {number} lod - Length Of Day in seconds
   */
  lod = 37;

  /**
   * @public
   * @member {number} xp - Pole x coordinate in arcseconds, call setEop() if member is changed without the constructor
   */
  xp = 0.0;
  /**
   * @public
   * @member {number} yp - Pole y coordinate in arcseconds, call setEop() if member is changed without the constructor
   */
  yp = 0.0;

  /**
   * @public
   * @member {number} dat - TAI - UTC in seconds, call setEop() if member is changed without the constructor
   */
  dat = 0.0;

  /**
   * @public
   * @member {number} dut1 - UT1 - UTC in seconds, call setEop() if member is changed without the constructor
   */
  dut1 = 0.0;

  /**
   * @public
   * @member {Vector} recef - Position in ECEF (ITRF) in km, call ecef2eci(t) to update eci counterpart
   */
  recef = vlib.vec(0, 0, 0);
  /**
   * @public
   * @member {Vector} vecef - Velocity in ECEF (ITRF) in km/s, call ecef2eci(t) to update eci counterpart
   */
  vecef = vlib.vec(0, 0, 0);
  /**
   * @public
   * @member {Vector} aecef - Acceleration in ECEF (ITRF) in km/s^2, call ecef2eci(t) to update eci counterpart
   */
  aecef = vlib.vec(0, 0, 0);

  /**
   * @public
   * @member {Vector} reci - Position in ECI (GCRS) in km, call eci2ecef(t) to update ecef counterpart
   */
  reci = vlib.vec(0, 0, 0);
  /**
   * @public
   * @member {Vector} veci - Velocity in ECI (GCRS) in km/s, call eci2ecef(t) to update ecef counterpart
   */
  veci = vlib.vec(0, 0, 0);
  /**
   * @public
   * @member {Vector} reci - Acceleration in ECI (GCRS) in km/s^2, call eci2ecef(t) to update ecef counterpart
   */
  aeci = vlib.vec(0, 0, 0);
  /**
   * @public
   * @member {Matrix3x3} Aeci2ecef - ECEF (ITRF) to ECI (GCRS) transformation matrix, call ecef2eci(t) or eci2ecef(t) to update this member, do not set value manually
   */
  Aecef2eci = vlib.eye(3);
  /**
   * @public
   * @member {Matrix3x3} Aeci2ecef - ECI (GCRS) to ECEF (ITRF) transformation matrix, call ecef2eci(t) or eci2ecef(t) to update this member, do not set value manually
   */
  Aeci2ecef = vlib.eye(3);

  // private members:
  #ass0;
  #a0si;
  #apni;
  #apn;
  #agst;
  #agsti;
  #ttt; // julian century of terrestrial time
  #jdut1; // julian day of ut1
  #xp; // in radians
  #yp; // in radians
  #omegaearth = vlib.vec(0, 0, 0);
  #t; // unix time in second UTC0
  #Aecef2pef;
  #Apef2ecef;
  #Apef2eci;
  #Aeci2pef;
  /**
   * Initialize Earth Orientation parameters at given time,
   * if time is not present, initialize EOP with default values
   * @param {Date} unixMs - Date value or unix time milliseconds in UTC0.
   * @param {String} relPath - Relative path of Earth Orientation folder w.r.t. project
   */
  constructor(unixMs, relPath) {
    if (relPath == undefined)
      folderPath = "./";
    else
      folderPath = relPath;
    this.#initialize_constants();
    this.setEopTime(unixMs);
  }

   /**
    * Return relative path w.r.t caller path, the path can only be set on new object declaration
   * @public
   * @member getRelativePath - relative path of EarthOrientation folder w.r.t caller path, only used during class object construct, do not set
   * @returns - relative folder path to .dat and .csv files
   */
   getRelativePath(){ return folderPath };
   
  /**
   * Update Earth Orientation Parameters (EOP) of the class manually
   * If this method called without argument, it uses the public EOP values to initialize private members
   * @public
   * @param {Object} eop -  EOP object
   * @param {number} eop.xp -  Pole x coordinate in arcseconds, call setEop() if member is changed without the constructor
   * @param {number} eop.yp -  Pole y coordinate in arcseconds, call setEop() if member is changed without the constructor
   * @param {number} eop.dut1 -  UT1 - UTC in seconds, call setEop() if member is changed without the constructor
   * @param {number} eop.lod -  Length Of Day in seconds
   * @param {number} eop.dat -  TAI - UTC in seconds, call setEop() if member is changed without the constructor
   */
  setEop(eop) {
    if (eop != undefined) {
      this.xp = eop.xp;
      this.yp = eop.yp;
      this.dut1 = eop.dut1;
      this.lod = eop.lod;
      this.dat = eop.dat;
    }
    const conv = Math.PI / (180 * 3600); // arcseconds to radians
    this.#xp = this.xp * conv;
    this.#yp = this.yp * conv;

    // perform time conversions
    const jd = unix2jd(this.#t);
    this.#jdut1 = jd + this.dut1 / 86400.0;
    // const jdtai = jd + dat/86400.0; // jd of atomic time
    const jdtt = jd + (this.dat + 32.184) / 86400.0; // jd of terrestrial time
    this.#ttt = (jdtt - 2451545.0) / 36525.0;
  }

  /**
   * Updates Earth Orientation Parameters (EOP) of given time, using EOP data table retrieved from
   * https://celestrak.org/SpaceData/EOP-Last5Years.csv if the stored table is out of date.
   * Then updates Aecef2eci and Aeci2ecef properties upon call
   * needs to be called when precise transformation is necessary (e.g., when datetime of subsequent calls to eci2ecef or ecef2eci are more than one week)
   * NOTE: the constructor calls this function initially, therefore within close time to the constructor date time,  user may not need to call this method
   * @public
   * @param {tUnixMs} eop Date value or unix time milliseconds in UTC0.
   */
  setEopTime(tUnixMs) {
    this.#t = tUnixMs / 1000.0;
    checkAndUpdateEopFile(this.#t);
    this.setEop(readEopParameters(this.#t));
    this.transform(tUnixMs);
  }

  // initialize constant variables #ass0; #a0si; #apni; #apn; #agst; #agsti
  // from corresponding .dat files
  #initialize_constants() {
    // % " to rad
    const convrtu = (0.000001 * Math.PI) / (180.0 * 3600.0); // % if micro arcsecond
    const convrtm = (0.001 * Math.PI) / (180.0 * 3600.0); // % if milli arcsecond

    const as = file2Matrix(folderPath+"iau00s.dat");
    // multiply 2nd and 3rd columns with convrtu and assign to ass0
    const ass0 = as
      .map((row) => row.slice(1, 3))
      .map((row) => row.map((col) => col * convrtu));
    const a0si = as.map((row) => row.slice(3));

    const an = file2Matrix(folderPath+"iau03n.dat");
    // first 5 columns are assigned to apni
    const apni = an.map((row) => row.slice(0, 5));
    // rest of collumns are assigned to apn
    const apn = an
      .map((row) => row.slice(6))
      .map((row) => row.map((col) => col * convrtm));

    // read agst and agsti from iau00gs.dat file:
    const ag = file2Matrix(folderPath+"iau00gs.dat");
    const agst = ag
      .map((row) => row.slice(1, 3))
      .map((row) => row.map((col) => col * convrtu));
    const agsti = ag.map((row) => row.slice(3));
    // fs.writeFileSync("pnConstants2.txt", JSON.stringify({agsti:agsti, agst:agst}), "utf8");

    this.#a0si = a0si;
    this.#ass0 = ass0;
    this.#apn = apn;
    this.#apni = apni;
    this.#agst = agst;
    this.#agsti = agsti;
  }

  /**
   * TODO: to be implemented
   * @param {Date} unixMs Date value or unix time milliseconds in UTC0.
   * @returns {EarthOrientation}
   */
  ecef2eci(unixMs) {
    this.transform(unixMs);
  }

  
  /**
   * calculates Aecef2eci and Aeci2ecef transformation matrix, does not update ECEF or ECI pos, vel and acc values
   * Implemented using Vallado's iau00f2i.m file source code with IAU2000b model
   * 
   * @param {Date} unixMs Date value or unix time milliseconds in UTC0.
   * @returns {EarthOrientation} the EarthOrientation object
   */
  transform(unixMs) {
    this.#t = unixMs / 1000.0;
    const ttt = this.#ttt;
    // %                           function iau00pnb
    // %
    // %  this function calulates the transformation matrix that accounts for the
    // %  effects of precession-nutation in the iau2000b theory.
    // %
    // %  reference    :
    // %    vallado       2004, 212-214

    // % " to rad
    const convrt = Math.PI / (180.0 * 3600.0);

    // % obtain data for calculations form the 2000b theory
    const args = fundarg(ttt);
    let pnsum = 0.0;
    let ensum = 0.0;
    for (let i = 76; i >= 0; i--) {
      const tempval =
        this.#apni[i][0] * args.l +
        this.#apni[i][1] * args.l1 +
        this.#apni[i][2] * args.f +
        this.#apni[i][3] * args.d +
        this.#apni[i][4] * args.omega;
      pnsum +=
        (this.#apn[i][0] + this.#apn[i][1] * ttt) * Math.sin(tempval) +
        (this.#apn[i][4] + this.#apn[i][5] * ttt) * Math.cos(tempval);
      ensum +=
        (this.#apn[i][2] + this.#apn[i][3] * ttt) * Math.cos(tempval) +
        (this.#apn[i][6] + this.#apn[i][7] * ttt) * Math.sin(tempval);
    }

    // % ------ form the planetary arguments
    let pplnsum = -0.000135 * convrt; //% " to rad
    let eplnsum = 0.000388 * convrt;

    // %  add planetary and luni-solar components.
    let deltapsi = pnsum + pplnsum;
    let deltaeps = ensum + eplnsum;
    const p = precess(ttt);

    const oblo = 84381.406 * convrt; //% " to rad

    // % ----------------- find nutation matrix ----------------------
    // % mean to true
    const a1 = vlib.rot1mat(p.ea + deltaeps);
    const a2 = vlib.rot3mat(deltapsi);
    const a3 = vlib.rot1mat(-p.ea);

    // % j2000 to date (precession)
    const a4 = vlib.rot3mat(-p.xa);
    const a5 = vlib.rot1mat(p.wa);
    const a6 = vlib.rot3mat(p.psia);
    const a7 = vlib.rot1mat(-oblo);

    // % icrs to j2000
    const a8 = vlib.rot1mat(-0.0068192 * convrt);
    const a9 = vlib.rot2mat(0.041775 * Math.sin(oblo) * convrt);
    // ?%      a9  = rot2mat(0.0166170*convrt);
    const a10 = vlib.rot3mat(0.0146 * convrt);

    // a10*a9*a8*a7*a6*a5*a4
    const prec = vlib.MxM(
      a10,
      vlib.MxM(a9, vlib.MxM(a8, vlib.MxM(a7, vlib.MxM(a6, vlib.MxM(a5, a4)))))
    );

    // a3*a2*a1
    const nut = vlib.MxM(a3, vlib.MxM(a2, a1));

    const pnb = vlib.MxM(prec, nut);
    // %                     END OF function iau00pnb

    //************************************************************************** */

    // %                           function iau00gst
    // %
    // %  this function finds the iau2000 greenwich sidereal time.
    // %
    // %  outputs       :
    // %    gst         - greenwich sidereal time        0 to twopi rad
    // %    st          - transformation matrix
    // %
    // %  references    :
    // %    vallado       2004, 216

    const deg2rad = Math.PI / 180.0;

    const ttt2 = ttt * ttt;
    const ttt3 = ttt2 * ttt;
    const ttt4 = ttt2 * ttt2;

    // % mean obliquity of the ecliptic
    let epsa = 84381.448 - 46.84024 * ttt - 0.00059 * ttt2 + 0.001813 * ttt3; // % "
    epsa = (epsa / 3600.0) % 360.0; // % deg
    epsa = epsa * deg2rad; // % rad

    // %  evaluate the ee complementary terms
    let gstsum0 = 0.0;
    for (let i = 32; i >= 0; i--) {
      const tempval =
        this.#agsti[i][0] * args.l +
        this.#agsti[i][1] * args.l1 +
        this.#agsti[i][2] * args.f +
        this.#agsti[i][3] * args.d +
        this.#agsti[i][4] * args.omega +
        this.#agsti[i][5] * args.lonmer +
        this.#agsti[i][6] * args.lonven +
        this.#agsti[i][7] * args.lonear +
        this.#agsti[i][8] * args.lonmar +
        this.#agsti[i][9] * args.lonjup +
        this.#agsti[i][10] * args.lonsat +
        this.#agsti[i][11] * args.lonurn +
        this.#agsti[i][12] * args.lonnep +
        this.#agsti[i][13] * args.precrate;
      gstsum0 +=
        this.#agst[i][0] * Math.sin(tempval) +
        this.#agst[i][1] * Math.cos(tempval); // % rad
    }

    // gstsum1 = 0.0;
    // for j = 1: -1 : 1
    //     i = 33 + j;
    let i = 33;
    const tempval =
      this.#agsti[i][0] * args.l +
      this.#agsti[i][1] * args.l1 +
      this.#agsti[i][2] * args.f +
      this.#agsti[i][3] * args.d +
      this.#agsti[i][4] * args.omega +
      this.#agsti[i][5] * args.lonmer +
      this.#agsti[i][6] * args.lonven +
      this.#agsti[i][7] * args.lonear +
      this.#agsti[i][8] * args.lonmar +
      this.#agsti[i][9] * args.lonjup +
      this.#agsti[i][10] * args.lonsat +
      this.#agsti[i][11] * args.lonurn +
      this.#agsti[i][12] * args.lonnep +
      this.#agsti[i][13] * args.precrate;
    const gstsum1 =
      (this.#agst[i][0] * Math.sin(tempval) +
        this.#agst[i][1] * Math.cos(tempval)) *
      ttt;

    const eect2000 = gstsum0 + gstsum1 * ttt; // % rad

    // % equation of the equinoxes
    const ee2000 = deltapsi * Math.cos(epsa) + eect2000; //  % rad

    // %  earth rotation angle
    const tut1d = this.#jdut1 - 2451545.0;
    const twopi = 2 * Math.PI;
    let era = twopi * (0.779057273264 + 1.00273781191135448 * tut1d);
    era = era % twopi; //  % rad

    // %  greenwich mean sidereal time, iau 2000.
    const gmst2000 =
      era +
      (0.014506 +
        4612.15739966 * ttt +
        1.39667721 * ttt2 -
        0.00009344 * ttt3 +
        0.00001882 * ttt4) *
        convrt; // % " to rad

    const gst = gmst2000 + ee2000; // % rad

    // % transformation matrix
    const st = vlib.rot3mat(-gst);

    const pm = polarm(this.#xp, this.#yp, ttt);

    //% ---- setup parameters for velocity transformations
    const thetasa = 7.29211514670698e-5 * (1.0 - this.lod / 86400.0);

    // reci = Apef2eci * Aecef2pef*recef
    this.#omegaearth = vlib.vec(0, 0, thetasa);
    this.#Aecef2pef = pm;
    this.#Apef2ecef = vlib.MT(pm);
    this.#Apef2eci = vlib.MxM(pnb, st);
    this.#Aeci2pef = vlib.MT(this.#Apef2eci);
    this.Aecef2eci = vlib.MxM(this.#Apef2eci, this.#Aecef2pef);
    this.Aeci2ecef = vlib.MT(this.Aecef2eci);
    return this;
  }
}

// test code
// const t = new Date("2019-01-01T00:00:00.000Z");
// const EO = new EarthOrientation(t);
// console.log(EO.Aecef2eci);

module.exports = {
  EarthOrientation,
}
