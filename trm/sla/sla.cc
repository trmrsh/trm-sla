//
// Python/C interface file for sla-related stuff
//

#include <Python.h>
#include "numpy/arrayobject.h"
#include "slalib.h"
#include "slamac.h"
#include "trm_vec3.h"
#include "trm_constants.h"

#include <iostream>

// Implements slaDtt

static PyObject* 
sla_dtt(PyObject *self, PyObject *args)
{

    double utc;
    if(!PyArg_ParseTuple(args, "d:sla.dtt", &utc))
	return NULL;

    double d = slaDtt(utc);

    return Py_BuildValue("d", d);

};

// Implements slaCldj

static PyObject* 
sla_cldj(PyObject *self, PyObject *args)
{

    int year, month, day;
    if(!PyArg_ParseTuple(args, "iii:sla.cldj", &year, &month, &day))
	return NULL;

    double mjd;
    int status;
    slaCldj(year, month, day, &mjd, &status);
    if(status == 1){
	PyErr_SetString(PyExc_ValueError, ("sla.cldj: bad year = " + Subs::str(year)).c_str());
	return NULL;
    }else if(status == 2){
	PyErr_SetString(PyExc_ValueError, ("sla.cldj: bad month = " + Subs::str(month)).c_str());
	return NULL;
    }else if(status == 3){
	PyErr_SetString(PyExc_ValueError, ("sla.cldj: bad day = " + Subs::str(day)).c_str());
	return NULL;
    }
    return Py_BuildValue("d", mjd);

};


static PyObject* 
sla_djcl(PyObject *self, PyObject *args)
{

    double mjd;
    if(!PyArg_ParseTuple(args, "d", &mjd))
	return NULL;

    int year, month, day;
    double frac;
    int status;
    slaDjcl(mjd, &year, &month, &day, &frac, &status);
    if(status == -1){
	PyErr_SetString(PyExc_ValueError, ("sla.djcl: bad date, < 4701 BC March 1"));
	return NULL;
    }
    return Py_BuildValue("iiid", year, month, day, 24.*frac);
};


static PyObject* 
sla_eqgal(PyObject *self, PyObject *args)
{

    double ra, dec;
    if(!PyArg_ParseTuple(args, "dd", &ra, &dec))
	return NULL;

    // convert angles to those expected by sla routines
    const double CFAC = Constants::PI/180.;
    double rar    = CFAC*15.*ra;
    double decr   = CFAC*dec;
    
    double glong, glat;
    slaEqgal(rar, decr, &glong, &glat);
    glong /= CFAC;
    glat /= CFAC;
    return Py_BuildValue("dd", glong, glat);
};

static PyObject* 
sla_galeq(PyObject *self, PyObject *args)
{

    double l, b;
    if(!PyArg_ParseTuple(args, "dd", &l, &b))
	return NULL;

    // convert angles to those expected by sla routines
    const double CFAC = Constants::PI/180.;
    double lr   = CFAC*l;
    double br   = CFAC*b;
    
    double ra, dec;
    slaGaleq(lr, br, &ra, &dec);
    ra   /= 15.*CFAC;
    dec  /= CFAC;
    return Py_BuildValue("dd", ra, dec);
};

// Computes TDB time corrected for light travel given a UTC (in MJD), 
// a target position (ICRS), and observatory position

static PyObject* 
sla_utc2tdb(PyObject *self, PyObject *args)
{
    PyObject *iutc = NULL;
    double ra, dec, longitude, latitude, height;
    double pmra = 0., pmdec = 0., epoch = 2000., parallax = 0., rv = 0.;
    if(!PyArg_ParseTuple(args, "Oddddd|ddddd:sla.utc2tdb", 
			 &iutc, &longitude, &latitude, &height, &ra, &dec, 
			 &pmra, &pmdec, &epoch, &parallax, &rv))
	return NULL;

    // Some checks on the inputs
   // Some checks on the inputs
    npy_intp nutc = 0;
    double vutc = 0.;
    if(iutc){
        if(PyFloat_Check(iutc)){
            vutc = PyFloat_AsDouble(iutc);
            if(PyErr_Occurred()){
                PyErr_SetString(PyExc_ValueError, "sla.amass: could not translate utc value");
                return NULL;
            }
            nutc = 1;
        }else if(PyArray_Check(iutc)){
            int nd = PyArray_NDIM(iutc);
            if(nd != 1){
                PyErr_SetString(PyExc_ValueError, "sla.amass: utc must be a 1D array or a float");
                return NULL;
            }
        }else{
            PyErr_SetString(PyExc_TypeError, "sla.amass: utc must be a 1D array or a float");
            return NULL;
        }
    }else{
        PyErr_SetString(PyExc_ValueError, "sla.amass: utc not defined");
        return NULL;
    }

    if(longitude < -360. || longitude > +360.){
	PyErr_SetString(PyExc_ValueError, "sla.utc2tdb: longituge out of range -360 to +360");
	return NULL;
    }

    if(latitude < -90. || latitude > +90.){
	PyErr_SetString(PyExc_ValueError, "sla.utc2tdb: latitude out of range -90 to +90");
	return NULL;
    }

    if(ra < 0. || ra > 24.){
	PyErr_SetString(PyExc_ValueError, "sla.utc2tdb: ra out of range 0 to 24");
	return NULL;
    }

    if(dec < -90. || dec > +90.){
	PyErr_SetString(PyExc_ValueError, "sla.utc2tdb: declination out of range -90 to +90");
	return NULL;
    }

    // convert angles to those expected by sla routines
    const double CFAC = Constants::PI/180.;
    double latr   = CFAC*latitude;
    double longr  = CFAC*longitude;
    double rar    = CFAC*15.*ra;
    double decr   = CFAC*dec;
    double pmrar  = CFAC*pmra/3600.;
    double pmdecr = CFAC*pmdec/3600.;

    double u, v;
    slaGeoc( latr, height, &u, &v);
    u *= Constants::AU/1000.0;
    v *= Constants::AU/1000.0;
    if(nutc){
        double tt  = vutc + slaDtt(vutc)/Constants::DAY;
        double tdb = tt  + slaRcc(tt, vutc-int(vutc), -longr, u, v)/Constants::DAY;

        // Compute position of Earth relative to the centre of the Sun and
        // the barycentre
        double ph[3], pb[3], vh[3], vb[3];
        slaEpv(tdb, ph, vh, pb, vb);
        
        // Create 3 vectors for simplicity of code
        Subs::Vec3 hpos(ph), bpos(pb), hvel(vh), bvel(vb);
        
        // Calculate correction from centre of Earth to observatory
        double last = slaGmst(tdb) + longr + slaEqeqx(tdb);
        double pv[6];
        slaPvobs(latr, height, last, pv);
        
        // Correct for precession/nutation
        double rnpb[3][3];
        slaPneqx(tdb, rnpb);
        slaDimxv(rnpb, pv, pv);
        slaDimxv(rnpb, pv+3, pv+3);
        
        Subs::Vec3 padd(pv), vadd(pv+3);
        
        // heliocentric
        hpos += padd;
        hpos *= Constants::AU;
        hvel += Constants::DAY*vadd;
        hvel *= Constants::AU/Constants::DAY;

        // barycentric
        bpos += padd;
        bpos *= Constants::AU;
        bvel += Constants::DAY*vadd;
        bvel *= Constants::AU/Constants::DAY;
        
        // At this point 'hpos' and 'bpos' contains the position of the observatory on the 
        // BCRS reference frame in metres relative to the helio- and barycentres. Now update 
        // the target position using space motion data.
        double nepoch = slaEpj(vutc);
        slaPm(rar, decr, pmrar, pmdecr, parallax, rv, epoch, nepoch, &rar, &decr);
        
        // Compute position vector of target
        double tv[3];
        slaDcs2c(rar, decr, tv);
        Subs::Vec3 targ(tv);
        
        // Finally, the helio- and barycentrically corrected times 
        double hcorr = dot(targ, hpos)/Constants::C/Constants::DAY;
        double bcorr = dot(targ, bpos)/Constants::C/Constants::DAY;
        
        double btdb  = tdb + bcorr;
        double htdb  = tdb + hcorr;
        double hutc  = vutc + hcorr;
        
        // and the radial velocities
        double vhel = -dot(targ, hvel)/1000.;
        double vbar = -dot(targ, bvel)/1000.;
        return Py_BuildValue("ddddddd", tt, tdb, btdb, hutc, htdb, vhel, vbar);

    }else{

        // an array has been passed and we need to do lots of setting up
        // An array has been passed, lots of painful setup ...
        nutc = PyArray_Size(iutc);
        PyObject *arr = PyArray_FROM_OTF(iutc, NPY_DOUBLE, NPY_IN_ARRAY);
        if(arr == NULL) return NULL;
        double *utc = (double *)PyArray_DATA(arr);

        npy_intp dim[1] = {nutc};
        PyArrayObject *ott = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
        if(ott == NULL){
            Py_DECREF(arr);
            return NULL;
        }
        PyArrayObject *otdb = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
        if(otdb == NULL){
            Py_DECREF(arr);
            Py_DECREF(ott);
            return NULL;
        }
        PyArrayObject *obtdb = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
        if(obtdb == NULL){
            Py_DECREF(arr);
            Py_DECREF(ott);
            Py_DECREF(otdb);
            return NULL;
        }
        PyArrayObject *ohutc = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
        if(ohutc == NULL){
            Py_DECREF(arr);
            Py_DECREF(ott);
            Py_DECREF(otdb);
            Py_DECREF(obtdb);
            return NULL;
        }
        PyArrayObject *ohtdb = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
        if(ohtdb == NULL){
            Py_DECREF(arr);
            Py_DECREF(ott);
            Py_DECREF(otdb);
            Py_DECREF(obtdb);
            Py_DECREF(ohutc);
            return NULL;
        }
        PyArrayObject *ovhel = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
        if(ovhel == NULL){
            Py_DECREF(arr);
            Py_DECREF(ott);
            Py_DECREF(otdb);
            Py_DECREF(obtdb);
            Py_DECREF(ohutc);
            Py_DECREF(ohtdb);
            return NULL;
        }
        PyArrayObject *ovbar = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
        if(ovbar == NULL){
            Py_DECREF(arr);
            Py_DECREF(ott);
            Py_DECREF(otdb);
            Py_DECREF(obtdb);
            Py_DECREF(ohutc);
            Py_DECREF(ohtdb);
            Py_DECREF(ovhel);
            return NULL;
        }
        
        // data pointers
        double *tt    = (double*) ott->data;
        double *tdb   = (double*) otdb->data;
        double *btdb  = (double*) obtdb->data;
        double *hutc  = (double*) ohutc->data;
        double *htdb  = (double*) ohtdb->data;
        double *vhel  = (double*) ovhel->data;
        double *vbar  = (double*) ovbar->data;
        
        for(int i=0; i<nutc; i++){
            tt[i]  = utc[i] + slaDtt(utc[i])/Constants::DAY;
            tdb[i] = tt[i]  + slaRcc(tt[i], utc[i]-int(utc[i]), -longr, u, v)/Constants::DAY;

            // Compute position of Earth relative to the centre of the Sun and
            // the barycentre
            double ph[3], pb[3], vh[3], vb[3];
            slaEpv(tdb[i], ph, vh, pb, vb);
        
            // Create 3 vectors for simplicity of code
            Subs::Vec3 hpos(ph), bpos(pb), hvel(vh), bvel(vb);
        
            // Calculate correction from centre of Earth to observatory
            double last = slaGmst(tdb[i]) + longr + slaEqeqx(tdb[i]);
            double pv[6];
            slaPvobs(latr, height, last, pv);
        
            // Correct for precession/nutation
            double rnpb[3][3];
            slaPneqx(tdb[i], rnpb);
            slaDimxv(rnpb, pv, pv);
            slaDimxv(rnpb, pv+3, pv+3);
        
            Subs::Vec3 padd(pv), vadd(pv+3);
        
            // heliocentric
            hpos += padd;
            hpos *= Constants::AU;
            hvel += Constants::DAY*vadd;
            hvel *= Constants::AU/Constants::DAY;

            // barycentric
            bpos += padd;
            bpos *= Constants::AU;
            bvel += Constants::DAY*vadd;
            bvel *= Constants::AU/Constants::DAY;
        
            // At this point 'hpos' and 'bpos' contains the position of the observatory on the 
            // BCRS reference frame in metres relative to the helio- and barycentres. Now update 
            // the target position using space motion data.
            double nepoch = slaEpj(utc[i]);
            slaPm(rar, decr, pmrar, pmdecr, parallax, rv, epoch, nepoch, &rar, &decr);
        
            // Compute position vector of target
            double tv[3];
            slaDcs2c(rar, decr, tv);
            Subs::Vec3 targ(tv);
        
            // Finally, the helio- and barycentrically corrected times 
            double hcorr = dot(targ, hpos)/Constants::C/Constants::DAY;
            double bcorr = dot(targ, bpos)/Constants::C/Constants::DAY;
        
            btdb[i]  = tdb[i] + bcorr;
            htdb[i]  = tdb[i] + hcorr;
            hutc[i]  = utc[i] + hcorr;
        
            // and the radial velocities
            vhel[i] = -dot(targ, hvel)/1000.;
            vbar[i] = -dot(targ, bvel)/1000.;

        }

        Py_DECREF(arr);

        return Py_BuildValue("OOOOOOO", ott, otdb, obtdb, ohutc, ohtdb, ovhel, ovbar);
    }

};

// Computes observational parameters such as airmass, altititude and elevation

static PyObject* 
sla_amass(PyObject *self, PyObject *args)
{

    PyObject *iutc = NULL;
    double ra, dec, longitude, latitude, height;
    double wave=0.55, pmra = 0., pmdec = 0., epoch = 2000., parallax = 0., rv = 0.;
    if(!PyArg_ParseTuple(args, "Oddddd|dddddd:sla.amass", 
			 &iutc, &longitude, &latitude, &height, &ra, &dec, 
			 &wave, &pmra, &pmdec, &epoch, &parallax, &rv))
	return NULL;

    // Some checks on the inputs
    npy_intp nutc = 0;
    double vutc = 0.;
    if(iutc){
        if(PyFloat_Check(iutc)){
            vutc = PyFloat_AsDouble(iutc);
            if(PyErr_Occurred()){
                PyErr_SetString(PyExc_ValueError, "sla.amass: could not translate utc value");
                return NULL;
            }
            nutc = 1;
        }else if(PyArray_Check(iutc)){
            int nd = PyArray_NDIM(iutc);
            if(nd != 1){
                PyErr_SetString(PyExc_ValueError, "sla.amass: utc must be a 1D array or a float");
                return NULL;
            }
        }else{
            PyErr_SetString(PyExc_TypeError, "sla.amass: utc must be a 1D array or a float");
            return NULL;
        }
    }else{
        PyErr_SetString(PyExc_ValueError, "sla.amass: utc not defined");
        return NULL;
    }

    if(longitude < -360. || longitude > +360.){
	PyErr_SetString(PyExc_ValueError, "sla.amass: longituge out of range -360 to +360");
	return NULL;
    }

    if(latitude < -90. || latitude > +90.){
	PyErr_SetString(PyExc_ValueError, "sla.amass: latitude out of range -90 to +90");
	return NULL;
    }

    if(ra < 0. || ra > 24.){
	PyErr_SetString(PyExc_ValueError, "sla.amass: ra out of range 0 to 24");
	return NULL;
    }

    if(dec < -90. || dec > +90.){
	PyErr_SetString(PyExc_ValueError, "sla.amass: declination out of range -90 to +90");
	return NULL;
    }

    if(wave <= 0. || wave > 1000000.){
	PyErr_SetString(PyExc_ValueError, "sla.amass: declination out of range -90 to +90");
	return NULL;
    }

    // convert angles to those expected by sla routines
    const double CFAC = Constants::PI/180.;
    double latr   = CFAC*latitude;
    double longr  = CFAC*longitude;
    double rar    = CFAC*15.*ra;
    double decr   = CFAC*dec;
    double pmrar  = CFAC*pmra/3600.;
    double pmdecr = CFAC*pmdec/3600.;

    // first three small corrections factors are assumed zero for ease of use
    const double DUT  = 0.; // UT1-UTC, seconds
    const double XP   = 0.; // polar motion, radians
    const double YP   = 0.; // polar motion, radians
    const double T    = 285.; // ambient temperature K
    const double P    = 1013.25; // ambient pressure, mbar
    const double RH   = 0.2;     // relative humidity (0-1)
    const double TLR  = 0.0065;  // lapse rate, K/metre

    double nepoch, zdob, decob, raob, refa, refb, tanz;
     
    if(nutc){
        // A single utc has been passed and we will return floats.
        double vazob, vhaob, vdelz, valtob, vairmass, vpaob;

        // correct for space motion
        nepoch = slaEpj(vutc);
        slaPm(rar, decr, pmrar, pmdecr, parallax, rv, epoch, nepoch, &rar, &decr);
            
        // *observed* azimuth (N->E), zenith distance, hour angle, declination, ra (all radians)
        slaI2o(rar, decr, vutc, DUT, longr, latr, height, XP, YP, T, P, RH, wave, TLR, 
               &vazob, &zdob, &vhaob, &decob, &raob);
        
        // compute refraction
        slaRefcoq(T, P, RH, wave, &refa, &refb);
        tanz = tan(zdob);
        vdelz = tanz*(refa + refb*tanz*tanz)/CFAC;
            
        // convert units
        valtob   = 90.-zdob/CFAC;
        vairmass = slaAirmas(zdob); 
        vazob   /= CFAC;
            
        // Compute pa
        vpaob = slaPa(vhaob,decr,latr)/CFAC;
        vpaob = vpaob > 0. ? vpaob : 360.+vpaob;
            
        vhaob *= 24./Constants::TWOPI;

        // return  airmass, altitude (deg), azimuth (deg, N=0, E=90),
        // hour angle, parallactic angle, angle of refraction
        return Py_BuildValue("dddddd", vairmass, valtob, vazob, vhaob, vpaob, vdelz);

    }else{

        // An array has been passed, lots of painful setup ...
        nutc = PyArray_Size(iutc);
        PyObject *arr = PyArray_FROM_OTF(iutc, NPY_DOUBLE, NPY_IN_ARRAY);
        if(arr == NULL) return NULL;
        double *utc = (double *)PyArray_DATA(arr);

        npy_intp dim[1] = {nutc};
        PyArrayObject *oair = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
        if(oair == NULL){
            Py_DECREF(arr);
            return NULL;
        }
        PyArrayObject *oalt = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
        if(oalt == NULL){
            Py_DECREF(arr);
            Py_DECREF(oair);
            return NULL;
        }
        PyArrayObject *oaz = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
        if(oaz == NULL){
            Py_DECREF(arr);
            Py_DECREF(oair);
            Py_DECREF(oalt);
            return NULL;
        }
        PyArrayObject *oha = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
        if(oha == NULL){
            Py_DECREF(arr);
            Py_DECREF(oair);
            Py_DECREF(oalt);
            Py_DECREF(oaz);
            return NULL;
        }
        PyArrayObject *opa = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
        if(opa == NULL){
            Py_DECREF(arr);
            Py_DECREF(oair);
            Py_DECREF(oalt);
            Py_DECREF(oaz);
            Py_DECREF(oha);
            return NULL;
        }
        PyArrayObject *odelz = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
        if(odelz == NULL){
            Py_DECREF(arr);
            Py_DECREF(oair);
            Py_DECREF(oalt);
            Py_DECREF(oaz);
            Py_DECREF(oha);
            Py_DECREF(opa);
            return NULL;
        }
        
        // data pointers
        double *airmass = (double*) oair->data;
        double *altob   = (double*) oalt->data;
        double *azob    = (double*) oaz->data;
        double *haob    = (double*) oha->data;
        double *paob    = (double*) opa->data;
        double *delz    = (double*) odelz->data;
        
        for(int i=0; i<nutc; i++){
            // correct for space motion
            nepoch = slaEpj(utc[i]);
            slaPm(rar, decr, pmrar, pmdecr, parallax, rv, epoch, nepoch, &rar, &decr);
            
            // *observed* azimuth (N->E), zenith distance, hour angle, declination, ra (all radians)
            slaI2o(rar, decr, utc[i], DUT, longr, latr, height, XP, YP, T, P, RH, wave, TLR, 
                   &azob[i], &zdob, &haob[i], &decob, &raob);
            
            // compute refraction
            slaRefcoq(T, P, RH, wave, &refa, &refb);
            tanz = tan(zdob);
            delz[i] = tanz*(refa + refb*tanz*tanz)/CFAC;
            
            // convert units
            altob[i]   = 90.-zdob/CFAC;
            airmass[i] = slaAirmas(zdob); 
            azob[i]   /= CFAC;
            
            // Compute pa
            paob[i] = slaPa(haob[i],decr,latr)/CFAC;
            paob[i] = paob[i] > 0. ? paob[i] : 360.+paob[i];
            
            haob[i] *= 24./Constants::TWOPI;
        }

        Py_DECREF(arr);

        // return  airmass, altitude (deg), azimuth (deg, N=0, E=90),
        // hour angle, parallactic angle, angle of refraction
        return Py_BuildValue("OOOOOO", oair, oalt, oaz, oha, opa, odelz);
    }
};


// Computes position of the Sun

static PyObject* 
sla_sun(PyObject *self, PyObject *args)
{

    double utc, longitude, latitude, height, wave=0.55, rh=0.2;
    int fast=1;
    if(!PyArg_ParseTuple(args, "dddd|ddi:sla.sun", &utc, &longitude, &latitude, &height, &wave, &rh, &fast))
	return NULL;

    // Some checks on the inputs
    if(longitude < -360. || longitude > +360.){
	PyErr_SetString(PyExc_ValueError, "sla.sun: longituge out of range -360 to +360");
	return NULL;
    }

    if(latitude < -90. || latitude > +90.){
	PyErr_SetString(PyExc_ValueError, "sla.sun: latitude out of range -90 to +90");
	return NULL;
    }

    if(wave <= 0. || wave > 1000000.){
	PyErr_SetString(PyExc_ValueError, "sla.sun: declination out of range -90 to +90");
	return NULL;
    }

    if(rh < 0. || rh > 1.){
	PyErr_SetString(PyExc_ValueError, "sla.sun: relative humidity out of range 0 to 1");
	return NULL;
    }

    // convert angles to those expected by sla routines
    const double CFAC = Constants::PI/180.;
    double latr   = CFAC*latitude;
    double longr  = CFAC*longitude;
    
    // first three small corrections factors are assumed zero for ease of use
    const double DUT  = 0.; // UT1-UTC, seconds
    const double T    = 285.; // ambient temperature K
    const double P    = 1013.25; // ambient pressure, mbar
    const double TLR  = 0.0065;  // lapse rate, K/metre

    // UT1
    double ut1 = utc + DUT;

    // TT
    double tt = utc + slaDtt(utc)/Constants::DAY;

    // Compute position of Earth relative to the centre of the Sun and
    // the barycentre
    double ph[3], pb[3], vh[3], vb[3];
    slaEvp(tt, -1., vb, pb, vh, ph);
    
    // Nutate
    double rmatn[3][3], phn[3];
    slaNut(tt, rmatn);
    slaDmxv(rmatn, ph, phn);

    // Calculate correction from centre of Earth to observatory
    double last = slaGmst(ut1) + longr + slaEqeqx(tt);
    double pv[6];
    slaPvobs(latr, height, last, pv);

    double sun[3];
    for(int i=0; i<3; i++)
	sun[i] = -phn[i]-pv[i];

    double ra, dec, az, el;
    slaDcc2s(sun, &ra, &dec);
    slaDe2h(last-ra, dec, latr, &az, &el);
    double refract=0.;

    if(fast){
	double refa, refb, zobs;
	slaRefcoq(T, P, rh, wave, &refa, &refb);
	slaRefz(Constants::PI/2.-el, refa, refb, &zobs);
	refract = Constants::PI/2.-el - zobs;
    }else{
	for(int i=0; i<5; i++){
	    double zd = Constants::PI/2.-el-refract;
	    slaRefro(zd, height, T, P, rh, wave, latr, TLR, 1.e-8, &refract);
	}
    }

    // return  azimuth and elevation
    return Py_BuildValue("ddddd", az/CFAC, (el+refract)/CFAC, refract/CFAC, 
			 12.*ra/Constants::PI, dec/CFAC);

};

// Convert FK4 B1950 to Fk5 J2000 coords

static PyObject* 
sla_fk425(PyObject *self, PyObject *args)
{

    double ra4, dec4, pmra4 = 0., pmdec4 = 0., parallax4 = 0., rv4 = 0.; 
    if(!PyArg_ParseTuple(args, "dd|dddd:sla.precess", &ra4, &dec4, &pmra4, 
			 &pmdec4, &parallax4, &rv4))
	return NULL;

    // Some checks on the inputs
    if(ra4 < 0. || ra4 > 24.){
	PyErr_SetString(PyExc_ValueError, "sla.fk425: ra out of range 0 to 24");
	return NULL;
    }

    if(dec4 < -90. || dec4 > +90.){
	PyErr_SetString(PyExc_ValueError, "sla.fk425: declination out of range -90 to +90");
	return NULL;
    }

    // convert angles to those expected by sla routines
    const double CFAC = Constants::PI/180.;
    double rar4    = CFAC*15.*ra4;
    double decr4   = CFAC*dec4;
    double pmrar4  = CFAC*pmra4/3600.;
    double pmdecr4 = CFAC*pmdec4/3600.;
    
    double rar5, decr5, pmrar5, pmdecr5, parallax5, rv5;
    slaFk425(rar4, decr4, pmrar4, pmdecr4, parallax4, rv4, 
	     &rar5, &decr5, &pmrar5, &pmdecr5, &parallax5, &rv5);

    double ra5    = rar5/CFAC/15.;
    double dec5   = decr5/CFAC;
    double pmra5  = 3600.*pmrar5/CFAC;
    double pmdec5 = 3600.*pmdecr5/CFAC;

    return Py_BuildValue("dddddd", ra5, dec5, pmra5, pmdec5, parallax5, rv5);

};

//----------------------------------------------------------------------------------------
// The methods

static PyMethodDef SlaMethods[] = {

    {"dtt", sla_dtt, METH_VARARGS, 
     "d = dtt(utc) returns TT-UTC in seconds. UTC in MJD = JD-2400000.5."},

    {"cldj", sla_cldj, METH_VARARGS, 
     "mjd = cldj(year, month, day) returns the mjd of the Gregorian calendar date; MJD = JD-2400000.5."},

    {"djcl", sla_djcl, METH_VARARGS, 
     "year,month,day,hour = djcl(mjd) decomposes an mjd into more usual times."},

    {"eqgal", sla_eqgal, METH_VARARGS, 
     "glong,glat = eqgal(ra,dec) returns galactic coords (degress) given ra, dec (J2000) in hours and degrees."},

    {"fk425", sla_fk425, METH_VARARGS, 
     "ra,dec,pmra,pmdec,parallax,rv = fk425(ra,dec,pmra=0,pmdec=0,parallax=0,rv=0) converts FK4 B1950 to FK5 J2000\n."
     "ra and dec are in hours and degrees; proper motions are in arcsec/year (not seconds of RA);\n"
     "parallax is in arcsec and the radial velocity is in km/s.\n"},

    {"galeq", sla_galeq, METH_VARARGS, 
     "ra,dec = galeq(glong,glat) returns FK5 J2000 coords (hours,degrees) given galactic coords (degrees)."},

    {"utc2tdb", sla_utc2tdb, METH_VARARGS, 
     "tt,tdb,btdb,hutc,htdb,vhel,vbar =\n"
     "    utc2tdb(utc,longitude,latitude,height,ra,dec,pmra=0,pmdec=0,epoch=2000,parallax=0,rv=0).\n\n"
     "All times are in MJD. Longitude and latitude are in degrees, east positive; ra and dec are in\n"
     "hours and degrees; proper motions are in arcsec/year (not seconds of RA); parallax is in arcsec\n"
     "and the radial velocity is in km/s. tt is terrestrial time (once ephemeris time); tdb is\n"
     "barycentric dynamical time; btdb is the barycentric dynamical time corrected for light travel\n"
     "time, i.e. as observed at the barycentre of the Solar system, hutc is the utc corrected for light\n"
     "travel to the heliocentre (usual form); htdb is the TDB corrected for light travel to the\n"
     "heliocentre (unusual). vhel and vbar are the apparent radial velocity of the target in km/s\n"
     "owing to observer's motion in relative to the helio- and barycentres. If utc is a float then\n"
     "so too will the output. If it is an array, then so will the outputs be."},

    {"amass", sla_amass, METH_VARARGS, 
     "airmass, alt, az, ha, pa, delz =\n"
     "  amass(utc,longitude,latitude,height,ra,dec,wave=0.55,pmra=0,pmdec=0,epoch=2000,parallax=0,rv=0).\n\n"
     "utc is an MJD or an array of MJDs. Longitude and latitude are in degrees, east positive; ra and dec are in\n"
     "hours and degrees; the wavelength of observation wave is in microns; proper motions are in\n"
     "arcsec/year (not seconds of RA); parallax is in arcsec and the radial velocity is in km/s.\n\n"
     "airmass is the airmass; alt and az are the observed altitude and azimuth in degrees with azimuth\n"
     "measured North through East; ha is the observed hour angle in hours; pa is the position angle\n"
     "of a parallactic slit; delz is the angle of refraction in degrees."},

    {"sun", sla_sun, METH_VARARGS, 
     "azimuth,elevation,refract,ra,dec = sun(utc,longitude,latitude,height,wave=0.55,rh=0.2,fast=True).\n\n"
     "All times are in MJD. Longitude and latitude are in degrees, east positive;\n"
     "the wavelength of observation wave is in microns; rh is the relative humidity.\n\n"
     "azimuth is measured in degrees, North through East, elevation in degrees above the horizon.\n"
     "refract is the angle of refraction in degrees. ra and dec are the position of the Sun for\n"
     "the mean equator and equinox of the utc supplied, FK5. fast determines whether a fast or slow\n"
     "is used. The fast method is OK for >15 degrees above the horizon, but for accurate values\n"
     "below this you may want the slow method\n"},

    {NULL, NULL, 0, NULL} /* Sentinel */
};

PyMODINIT_FUNC
init_sla(void)
{
    (void) Py_InitModule("_sla", SlaMethods);
    import_array();
}
