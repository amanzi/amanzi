import numpy as np
import spline
import scipy.interpolate


class VanGenuchten(object):
    def __init__( self, m=None, alpha=0., sr=0.0, smoothing_interval=0.0, n=None, l=0.5, smoothing_interval_sat=0.0):
        self._alpha = alpha
        self._sr = sr
        self._l = l

        if n is None:
            self._m = m
            self._n = 1.0/(1.0-m)
        else:
            self._m = 1 - 1.0/n 
            self._n = n

            
        self.s_r = self._sr

        # smoothing for kr
        self._s0 = 1.0 - smoothing_interval
        if self._s0 < 1.:
            self._spline = spline.Spline(self._s0, self.k_relative(self._s0), self.d_k_relative(self._s0),
                                         1.0, 1.0, 0.)
            
        # smoothing for sat
        self._pc0 = smoothing_interval_sat
        if self._pc0 > 0.:
            self._spline_sat = spline.Spline(0., 1., 0., self._pc0, self.saturation(self._pc0), self.d_saturation(self._pc0))

    def capillaryPressure( self, s ):
        if s <= self._sr:
            return np.inf
        if s >= 1.:
            return 0.
        se = (s - self._sr) / (1.0 - self._sr);
        if (se < 1.e-8):
            return pow(se, -1.0/(self._m * self._n)) / self._alpha
        else:
            return (pow(pow(se, -1.0/self._m) - 1.0, 1/self._n)) / self._alpha

    def d_capillaryPressure( self, s ):
        if s <= self._sr:
            return 0.
        if s >= 1.:
            return 0.
        se = (s - self._sr) / (1.0 - self._sr);
        if se < 1.e-8:
            return -1 / (self._m * self._n * self._alpha) \
              * pow(se, -1/(self._m * self._n) - 1.) / (1. - self._sr)
        else:
            return -1 / (self._m * self._n * self._alpha) \
              * pow(pow(se,-1./self._m) - 1., 1./self._n - 1.) \
              * pow(se, -1./self._m - 1.) / (1. - self._sr)

    def saturation( self, pc ):
        if pc <= 0.0:
            return 1.0
        elif pc < self._pc0:
            return self._spline_sat.Value(pc)
        else:
            se = pow(1.0 + pow(self._alpha*pc, self._n), -self._m)
            return  se * (1.0 - self._sr) + self._sr

    def d_saturation( self, pc ):
        if pc <= 0.:
            return 0.
        elif pc < self._pc0:
            return self._spline_sat.Derivative(pc)
        else:
            dse_dpc = -self._m*self._n * pow(1.0 + pow(self._alpha*pc, self._n), -self._m-1.0) * pow(self._alpha*pc, self._n-1) * self._alpha
            return  (1.0 - self._sr) * dse_dpc

    def k_relative( self, s ):
        if s >= 1.:
            return 1.
        elif s <= self._sr:
            return 0.
        elif s <= self._s0:
            se = (s - self._sr) / (1.0-self._sr)
            return (se**self._l) * pow( 1.0 - pow( 1.0 - pow(se,1.0/self._m),self._m), 2);
        else:
            return self._spline.Value(s)

    def d_k_relative( self, s ):
        if s >= 1.:
            return 0.
        elif s <= self._sr:
            return 0.
        elif s <= self._s0:
            se = (s - self._sr) / (1.0-self._sr)

            x = pow(se,1./self._m)
            y = pow(1. - x, self._m)
            dkdse = (1.0 - y) * ( (1. - y) / 2. + 2 * x * y / (1.-x)) / np.sqrt(se)
            return dkdse / (1 - self._sr)
        else:
            return self._spline.Derivative(s)

    def label( self ):
        return "VG: n=%1.2g,a=%1.2e,smooth=%g"%(self._n, self._alpha, 1-self._s0)
    
    
class VanGenuchtenSpline(object):
    def __init__( self, m=None, alpha=0., sr=0.0, smoothing_interval=0.0, n=None, l=0.5, n_control=100):
        self._vg = VanGenuchten(m,alpha,sr,smoothing_interval, n, l)

        #cs = np.linspace(0., 1., n_control)
        #
        cs = np.array([np.cos(i*np.pi/(n_control-1)) for i in range(n_control)]) 
        #cs = np.sign(cs) * np.abs(cs)**0.5  
        cs = (1 - cs)/2.

        cs = 1-cs
        cs = (cs * (1-sr) + sr)[:-1]
        
        pcs = np.array([self._vg.capillaryPressure(s) for s in cs])
        dsats = np.array([self._vg.d_saturation(pc) for pc in pcs])
       
        print "Control sats:", cs
        print "Control PCs:", pcs 
        print "Control dsats:", dsats
    
        krs = np.array([self._vg.k_relative(s) for s in cs])
        dkrs = np.array([self._vg.d_k_relative(s) for s in cs])
        
        print "Control krs:", krs
        print "Control dkrs:", dkrs
        
        # add beginning and end to avoid extrap, get derivs correct
        #self._control_sats = np.zeros((len(cs)+2,),'d'); self._control_sats[1:-1] = cs
        #self._control_pcs = np.zeros((len(cs)+2,),'d'); self._control_pcs[1:-1] = pcs
        #self._control_krs = np.zeros((len(cs)+2,),'d'); self._control_krs[1:-1] = krs
        
        #self._control_dsats = np.zeros((len(cs)+2,),'d'); self._control_dsats[1:-1] = dsats
        #self._control_dkrs = np.zeros((len(cs)+2,),'d'); self._control_dkrs[1:-1] = dkrs        

        #self._control_sats[0] = 1.1; self._control_sats[-1] = sr
        #self._control_dsats[0] = 0.; self._control_dsats[-1] = 0.
        #self._control_pcs[0] = -pcs[1]; self._control_pcs[-1] = 1.e12
        #self._control_krs[0] = 1.; self._control_krs[-1] = 0.
        #self._control_dkrs[0] = 0.; self._control_dkrs[-1] = 0.
            
        # splines
        #self._spline_sat = scipy.interpolate.PchipInterpolator(self._control_pcs, self._control_sats)
        #self._spline_kr = scipy.interpolate.PchipInterpolator(np.flipud(self._control_sats), np.flipud(self._control_krs))
        

        self._control_sats = cs
        self._control_pcs = pcs
        self._control_krs = krs
        self._control_dsats = dsats
        self._control_dkrs = dkrs
        
        self._spline_sat = spline.SplineList(self._control_pcs, self._control_sats, self._control_dsats, monotonicity_preserving=True)
        self._spline_kr = spline.SplineList(np.flipud(self._control_sats), np.flipud(self._control_krs), np.flipud(self._control_dkrs), monotonicity_preserving=True)

    
    def label(self):
        return "VGS:"+self._vg.label()[3:]
    
    def capillaryPressure( self, s ):
        raise RuntimeError("not implemented for spline")
    def d_capillaryPressure( self, s ):
        raise RuntimeError("not implemented for spline")
    
    def saturation( self, pc ):
        return self._spline_sat(pc)
    def d_saturation( self, pc ):
        #return (self._spline_sat(pc + 1) - self._spline_sat(pc-1))/2.
        return self._spline_sat.Derivative(pc)
    
    def k_relative( self, s ):
        return self._spline_kr(s)
    def d_k_relative( self, s ):
        #return (self._spline_kr(s+.0001) - self._spline_kr(s-0.0001))/0.0002
        return self._spline_kr.Derivative(s)

def xml_VanGenuchten(wrm, region):
    import sys, os
    from ats_python.utils import read_source_for_template
    from amanzi_xml.utils import search
    
    fname = os.path.join(os.environ['ATS_SRC_DIR'],
                         "src", "pks", "flow", "constitutive_relations", "wrm",
                         "wrm_van_genuchten.hh")

    template = read_source_for_template.readTemplate(fname)

    # fill the template
    template.attrib['name'] = region
    search.searchAndReplaceByName(template, "region=%s"%region)
    search.searchAndReplaceByName(template, "van Genuchten alpha=%g"%wrm._alpha)
    search.searchAndReplaceByName(template, "van Genuchten m=%g"%wrm._m)
    search.searchAndReplaceByName(template, "residual saturation=%g"%wrm._sr)
    search.searchAndReplaceByName(template, "smoothing interval width=%g"%wrm._pc0)
    return template
    
