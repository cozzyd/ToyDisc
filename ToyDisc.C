
/** Toy Disc Askaryan Neutrino Simulation 
 *  
 *  Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 *
 *  Implements sampling and a simple Askaryan radiation / detector model in uniform density block of ice. 
 *
 */

#include "TFile.h" 
#include "TTree.h" 
#include "TVector3.h" 
#include "TRandom3.h" 
#include "TMath.h" 
#include "TH2.h" 
#include "TF1.h" 
#include "TVectorD.h"

#include <fenv.h> 



const double freq_ref = 1.15e3;


/** this implements the cross-section model suggested in Connolly et al */ 
double cross_section(double E, bool cc, bool nubar) 
{
  
  const double C0[2][2] = { { -1.826, -1.826}, {-1.033, -1.033}}; 
  const double C1[2][2] = { { -17.31, -17.31}, {-15.95, -15.95}}; 
  const double C2[2][2] = { { -6.448, -6.406}, {-7.296, -7.247}}; 
  const double C3[2][2] = { { 1.431,1.431}, {1.569, 1.569}}; 
  const double C4[2][2] = { { -18.61,-17.91}, {-18.3, -17.72}}; 
   
  double e = log10(E)-9; 
  double lnterm = log(e-C0[nubar][cc]); 
  double log10_sigma = C1[nubar][cc] + 
                       C2[nubar][cc]*lnterm  +
                       C3[nubar][cc]*lnterm*lnterm +
                       C4[nubar][cc]/lnterm;

  return pow(10,log10_sigma); 
}


// ranodmly pick a flavor
int randomFlavor(TRandom * rng = gRandom) 
{
  int r = rng->Integer(6); 
  if (r < 3) return -1-r; 
  else return r-2; 
}




//randomly picks charged current or not
bool isCC(double E, TRandom * rng= gRandom) 
{
  double nc_frac = 0.252162 + 0.0256*log(log10(E)-9-1.76); //from Connolly et al. 
  return rng->Rndm() > nc_frac; 
}


//from icemc paper 
double coneWidth(bool em, double Eshower, double freq, double density = 0.917, double n = 1.78)
{

  if (em) 
  {
    const double Elpm = 2e15;
    return 2.7 * TMath::DegToRad() * freq_ref/freq *  pow(Elpm / (0.14 * Eshower + Elpm),0.3);
  }


  const double kdelta = 18.33; 
  const double x0 =40.3; 
  return 3e10/(freq*1e6) * density / kdelta / x0 /(n*n-1); 
}


//Connolly et al parametrization
double inelasticity(double E, bool cc, bool nubar, TRandom * rng = gRandom) 
{
  double e = log10(E)-9; 
  double R1 = rng->Rndm(); 

  bool low_y = R1 < (0.128) * sin(-0.197 * (e-21.8));

  const double A_low[4] = { 0, 0.0941, 4.72, 0.456};

  // cc, bar, index
  const double A_high[2][2][4] =
  { 
    {
      { -0.005, 0.23, 3, 1.7},
      { -0.005, 0.23, 3, 1.7}
    }, 
    {
      {-0.008, 0.26, 3, 1.7},
      {-0.0026, 0.085, 4.1, 1.7}
    }
  };

  const double * A = low_y ? A_low : A_high[cc][nubar]; 

  double C1 = A[0] - A[1] * exp(-(e-A[2])/A[3]);

  double R = rng->Rndm(); 
  if (low_y) 
  {
    double ymax = 1e-3;
    double ymin = 0; 
    double C2 = 2.55 -0.0949 * e; 
    return C1 + pow(R*pow(ymax-C1, -1./C2+1) + (1-R)* pow(ymin-C1,-1./C2+1), C2/(C2-1)); 
  }
  double ymax = 1; 
  double ymin = 1e-3; 

  return pow(ymax - C1,R) / pow(ymin-C1,R-1) + C1; 
}


TVector3 random_dir(double cos_theta_min = -1, double cos_theta_max = 1, double phi_min = 0, double phi_max = 2*TMath::Pi(), TRandom * rng = gRandom)
{
  double cos_theta = rng->Uniform(cos_theta_min, cos_theta_max); 
  if (phi_max < phi_min) phi_max += 2 *TMath::Pi(); 
  double phi = rng->Uniform(phi_min, phi_max); 
  double sin_theta = sqrt(1-cos_theta*cos_theta); 
  return TVector3( cos(phi) * sin_theta, sin(phi) *sin_theta, cos_theta); 
}

struct ToyDiscEvent
{
  TVector3 vertex;
  TVector3 direction;
  int flavor = 0; 
  double E = 0; 
  double Eshower = 0;
  double absorption_weight = 0; 
  double projected_area = 0; 
  double interaction_weight = 1; 
  double lint = 0;
  double Efield = 0; 
  double Ech_1m = 0; 
  double particle_depth = 0;
  double view_angle = 0; 
  double entrance_distance = 0; 
  double atten_factor = 0; 
  int nint =0; 
  bool pass = false ; 
  bool should_save = false; 
  bool forced_interaction = true; 
  bool in_volume = false; 
  int cc = 0; 
  double lproj = 0; 
}; 

// random exponential variate 
double ranexp(double tau, double min_val = 0, double max_val = -1, TRandom * rng = gRandom) 
{
 return -tau * log(rng->Uniform(max_val > 0 ? exp(-max_val/tau) : 0,exp(-min_val/tau))) ;
}


struct ToyDiscSetup
{
  double R = 15e3; 
  double h = 3e3; 
  double rho = 0.917; //g/cm^3
  double n = 1.78; 
  double freq_min = 150; //frequency used for askaryan calculation
  double freq_max = 600; //frequency used for askaryan calculation
  double freq_step = 10; 
  double threshold=10e-6; //V/m
  double save_thresold=1e-6; //V/m // save events with at least this much efield 
  double Rearth = 6371e3; 
  double threshold_rms = 1e-9; 
  double attenuation_length = 1e3; //m 
  int nbins_slant = 360; 
  int nbins_depth = 100; 
  TVector3 detector = TVector3(0,0,-100); //detector location
};

ToyDiscSetup default_setup; 


double prem_slant_depth(const TVector3 & pos, const TVector3 & dir, double ice_height = 3000, double earth_radius = 6371e3, double ice_density = 0.917); 

void setup_h_slant(TH2 & hslant, const ToyDiscSetup & s) 
{
  

  //let's create the histogram 

  for (int j = 1; j <= hslant.GetNbinsY(); j++) 
  {

    TVector3 pos(0,0,hslant.GetYaxis()->GetBinCenter(j)); 
    for (int i = 1; i <= hslant.GetNbinsX(); i++) 
    {
      double ang = hslant.GetXaxis()->GetBinCenter(i); 
      TVector3 dir(0, sin(ang), cos(ang)); 
      hslant.SetBinContent(i,j, prem_slant_depth(pos,dir,s.h, s.Rearth, s.rho)); 
    }

    hslant.SetTitle("Slant depth"); 
    hslant.GetXaxis()->SetTitle("Angle below horizontal [rad]"); 
    hslant.GetYaxis()->SetTitle("Depth [m]"); 
    hslant.GetZaxis()->SetTitle("Slant Depth [g/cm^{2}]"); 
    hslant.SetStats(0); 
  }

}


class ToyDiscSim  
{

  ToyDiscSetup s; 
  TRandom * rng ;
  double R2,h2; 

  double the_sampling_area;
  double the_sampling_volume; 
  double the_sampling_half_length; 
  TH2D hslant; 

  public: 
    ToyDiscSim(const ToyDiscSetup & setup = default_setup,  TRandom * rng = gRandom, int nbins_slant = 360) 
      : s(setup), rng(rng), hslant("hslant","SlantDepth;angle;grammage (g cm^2)", s.nbins_slant,0,TMath::Pi(), 1+s.nbins_depth, -s.h+(-s.h/s.nbins_depth)/2, s.h/s.nbins_depth/2) 
    {

      setup_h_slant(hslant,s); 

      R2 = s.R*s.R;
      h2 = s.h*s.h;

      the_sampling_area = 4* R2  + h2; 
      the_sampling_half_length = sqrt(R2+h2/4); 
      the_sampling_volume = TMath::Pi() * R2* s.h; 

    }


    const ToyDiscSetup * setup() { return &s; } 
    double sampling_area() const { return the_sampling_area; }
    double sampling_box_half_length() { return the_sampling_half_length ;} 
    double sampling_volume() const { return the_sampling_volume; } 


    
    void randomEvent(ToyDiscEvent* event, double E, bool forced = true, 
                     double cos_theta_min = -1, double cos_theta_max =1, double phi_min = 0, double phi_max = 2*TMath::Pi())
    {
      TVector3 dir = random_dir(cos_theta_min, cos_theta_max, phi_min, phi_max, rng); 
      TVector3 x = forced ? randomVertex() : randomTrajectory(dir); 
      doEvent(event,E, x,dir,forced); 
    }


    const TH2D & getHSlant() { return hslant ; }

    //returns true if in the volume! 
    bool doEvent(ToyDiscEvent * event, double E, const TVector3 & x0, const TVector3 & direction, bool forced_interaction = true, int flavor = 0, int cc = -1)
    {

      TVector3 x1, x2; 
      double ts[2]; 
      int nint = intersections(x0, direction, x1,x2,ts); 


      double theta = direction.Theta(); 
      event->projected_area = TMath::Pi() * R2 *fabs(cos(theta)) + 2* s.R * s.h * fabs(sin(theta)); 

      event->nint = nint; 
      event->E = E; 
      event->forced_interaction = forced_interaction;

      event->direction = direction; 
      if (nint == 0 || (forced_interaction && !inside(x0)))  //this misses the volume 
      {
        event->in_volume = false; 
        event->direction = direction; 
        event->vertex =  x0; 
        event->cc = -1; 
        event->flavor = 0; 
        event->absorption_weight = 0;
        event->interaction_weight  = 0; 
        event->particle_depth = 0;
        event->entrance_distance  = 0; 
        event->pass = false; 
        event->should_save = false; 
        event->Efield = 0;
        event->lint = 0; 
        event->lproj = 0; 
        event->view_angle = 0; 
        event->Ech_1m = 0; 
        event->atten_factor = 0; 
        event->Eshower = 0; 
        return false;
      }


      
      if (cc < 0) cc = isCC(E, rng); 
      if (!flavor) flavor = randomFlavor(rng); 

      TVector3 vertex = x0; 

      event->lproj = fabs(ts[1]-ts[0]); 
      double sigma = cross_section(E,cc,flavor < 0); 
      double number_density = s.rho / 1.66e-24; 
      event->lint = 1./(number_density*sigma*100); //m 

      if (!forced_interaction) 
      {
        double where = ranexp(event->lint, 0, event->lproj);
        vertex = (ts[0] < ts[1] ? x1 : x2)  + direction*where; 
        if (!inside(vertex) || std::isnan(vertex.X()))
        {
          printf("%g %g %g\n", where, event->lint, event->lproj);
          vertex.Print(); 
        }
      }

      event->in_volume = true; 
      event->vertex = vertex; 
      event->flavor = flavor; 
      event->cc = cc; 
      event->interaction_weight = 1-exp(-event->lproj/event->lint); 


      auto x = s.detector-vertex; 
      double det_angle = direction.Angle(x); 
      double distance = x.Mag(); 
      double atten_factor = exp(-distance/s.attenuation_length) / distance; 
      bool em = cc && abs(flavor)==1; 
      double y  =inelasticity(E,cc,flavor < 0, rng); 
      double Eshow = em ? (1-y)*E : y*E; 
      event->Eshower = Eshow; 
      double cos_cher = 1/s.n; 
      double cher = acos(cos_cher); 
      double sin_cher = sin(cher); 
      double view_angle = fabs(cher-det_angle); 
      double E1m_cher = 0; 
      double Efield = 0; 
      for (double freq = s.freq_min; freq <= s.freq_max; freq+= s.freq_step )
      {
       double this_E1m_cher = s.freq_step * 2.53e-7 * Eshow/1e12 * freq/freq_ref  / (1 + pow(freq/freq_ref,1.44));
       E1m_cher +=  this_E1m_cher; 
       double width = coneWidth(em, Eshow, freq, s.rho, s.n); 
       Efield +=  this_E1m_cher * atten_factor  * sin(view_angle)/sin_cher * TMath::Gaus(view_angle,0,width); 
      }

      event->Efield = Efield; 
      event->Ech_1m = E1m_cher; 
      event->atten_factor = atten_factor; 
      event->view_angle = view_angle; 
      event->pass =  Efield > rng->Gaus(s.threshold,s.threshold_rms); 
      event->should_save =  Efield > s.save_thresold; 

      event->particle_depth = hslant.Interpolate(theta,vertex.Z()) /  1.66e-24; // convert g/cm^2 to nucleons / cm^2

      //subtract the particle depth within the interaction volume 
      TVector3 entrance = ts[0] < ts[1] ? x1 : x2; 
      event->entrance_distance = (entrance-vertex).Mag()*100;  //in cm
      event->particle_depth -= (number_density * event->entrance_distance); 
      if (event->particle_depth < 0) event->particle_depth = 0; 
      event->absorption_weight = exp(-event->particle_depth*sigma); 

    return true; 
  }



  int intersections(const TVector3 & x0, const TVector3 & p, TVector3 & int1, TVector3 & int2, double * ts=0) const
  {
   
    int nint = 0; 

    //check if it intersects the top plane (z = 0) or bottom plane (z = -h); 
    
    if (p.Z()!=0) 
    {
      double t= -x0.Z()/p.Z();

      double x = x0.X() + p.X()*t; 
      double y = x0.Y() + p.Y()*t; 

      if (x*x+y*y < R2)
      {
        nint++; 
        int1 = TVector3(x,y,0); 
        if (ts) ts[0] = t; 
      }

      t = (-s.h-x0.Z())/p.Z();

      x = x0.X() + p.X()*t; 
      y = x0.Y() + p.Y()*t; 
      if (x*x+y*y < R2)
      {
        nint++; 
        (nint == 1? int1 : int2) = TVector3(x,y,-s.h); 
        if (ts) ts[nint-1] = t; 
      }
    }


    //check if it intersects the sides 
    if (nint  < 2 ) 
    {


      double dx = p.X(); 
      double dy = p.Y(); 
      double dr2 = dx*dx+dy*dy; 

      double D = x0.X() * (x0.Y() +dy) - (x0.X()+dx)*x0.Y(); 


      double discr = R2*dr2 -D*D; 
      if (discr < 0)
      {
        return nint;
      }

      for (int sign = 1; sign >=-1; sign-=2) 
      {
        double xx = D * dy  + sign*(dy > 0 ? 1 : -1) * dx *sqrt(discr); 
        xx/=(dr2); 
        double yy = -D * dx  + sign*fabs(dy) *sqrt(discr); 
        yy/=(dr2);
        double t = (p.X() ==0) ?  (yy-x0.Y()/p.Y()) : (xx-x0.X())/p.X();

        double zz = x0.Z() + t * p.Z(); 

        if (zz < 0 && zz > -s.h) 
        {
          nint++; 
          (nint == 1 ? int1 : int2) = TVector3(xx,yy,zz); 
          if (ts) ts[nint-1] = t; 
        }

        if (discr == 0) break; 
      }
    }

    return nint; 
  }

  bool inside(const TVector3 & v)  const
  {
    if (v.Z() > 0) return false;
    if (v.Z() < -s.h) return false;
    if (v.Perp2() > R2) return false; 
    return true; 
  }


  TVector3 randomVertex()
  {
    double z = rng->Uniform(-s.h,0); 
    double x,y;
    do 
    {
      x = rng->Uniform(-s.R,s.R);
      y = rng->Uniform(-s.R,s.R);
    } while (x*x+y*y > R2); 

    return TVector3(x,y,z); 
  }

  // for a direction
  TVector3  randomTrajectory(const TVector3 & dir )
  {

    TVector3 u = dir.Orthogonal(); 
    TVector3 v=  dir.Cross(u); 
    TVector3 ref(0,0,-s.h/2); 
    double L = sampling_box_half_length(); 
    return ref + rng->Uniform(-L,L)*u + rng->Uniform(-L,L) * v; 
  }


};





//modified so that ice is near the top
double prem_density(double r, double ice_height = 3000, double earth_radius = 6371e3, double ice_density = 0.917) 
{

  if (r ==0) return  1.3088e1 ; 
  if (r >= earth_radius) return 0; 

  double rs[11] = { 1.2215e6, 3.48e6, 3.63e6, 5.701e6, 5.77e6, 5.971e6, 6.151e6, 6.3466e6, 6.356e6, earth_radius-ice_height, earth_radius};
  double as[11] = { -2.1773e-10, -2.4123e-10, 0,-3.0922e-11,0,0,0,0,0,0,0}; 
  double bs[11] = { 1.911e-8, 1.3976e-4, -5.0007e-4, -2.4441e-4, -2.3286e-4, -1.2603e-3, -5.9706e-4, 1.0869e-4,0,0,0};
  double cs[11] = { 1.3088e4, 1.2346e4, 7.3067e3, 6.7823e3, 5.3197e3, 1.1249e4, 7.1083e3, 2.691e3, 2.9e3, 2.6e3, ice_density*1e3};

  int i = std::upper_bound(rs, rs + 11, r) - rs; 

  return 1e-3 * ( as[i]*r*r + bs[i] *r + cs[i]); 
}



/* Integrate backwards from position until outside the surface of the earth. Returns slant depth in g/cm^2
 * */

TVector3 origin(0,0,0); 

int sphere_intersection(const TVector3 & pos, const TVector3 & dir, double * ts,  double R = 6371e3, const TVector3 & x0 = origin) 
{
  TVector3 o = pos-origin; 
  TVector3 l = dir.Unit(); 
  double lo = l*o; 
  double discr2 = lo*lo - (o.Mag2()-R*R); 
  int nint = discr2 < 0 ? 0 : 
            discr2 == 0 ? 1 : 
            2; 


  if (ts && nint) 
  {
    double discr = sqrt(discr2); 
    ts[0] = -lo - discr; 
    ts[1] = -lo + discr; 
  }

  return nint; 



}
  



double prem_slant_depth(const TVector3 & pos, const TVector3 & dir, 
                        double ice_height , double earth_radius , double ice_density ) 
{

  double grammage = 0; 
  TVector3 x0 = pos + TVector3(0,0,earth_radius); //move to top of ice 
  TVector3 v = dir.Unit(); 

  //find intersection with sphere 
  double ts[2]; 
  int nint = sphere_intersection(x0,dir, ts, earth_radius); 
  if (nint < 2) return 0; 

  double tmin = ts[0]; // this was the spehre intersection that the neturino saw before
  double tmax = TMath::Min(0.,ts[1]); // integrate either to us or earth surface, whichever was sooner

  //use TF1 interface for integral out of laziness
  TF1 f("fintegral", [=](double *xx, double *P) 
      {
      double t = *xx; 
      double r =  (x0 + t*v).Mag(); 
      return prem_density(r, ice_height, earth_radius, ice_density);
      }, tmin, tmax,0); 




  return f.Integral(tmin,tmax,1e-3)*100;  //100 to convert m to cm in integration
}



void ToyDisc(int forced = 1, double E = 1e18,  int N = 1e8, double cos_theta_min = -1,
              double cos_theta_max = 1, double phi_min = 0, double phi_max = TMath::Pi()*2 ) 
{

  printf("Running in %s mode at E=%g eV for %d events.\n", forced ? "FORCED" : "UNFORCED", E, N ); 
  double str = (cos_theta_max-cos_theta_min) * (phi_max-phi_min) / TMath::Pi(); 
  printf("Number of steradians is %g pi \n", str); 

  ToyDiscSim sim; 

  TFile f(Form("toydisc_%d_%g.root",forced, log10(E)),"RECREATE"); 

  TTree * t = new TTree("toy","toy"); 
  ToyDiscEvent * event = new ToyDiscEvent;
  t->Branch("event",&event); 
  int i ; 
  t->Branch("i",&i); 

  double sum_weights = 0; 
  int one_percent = N/100; 
  int last_percent= 0; 

  for (i = 0; i < N; i++) 
  {
    sim.randomEvent(event, E, forced, cos_theta_min, cos_theta_max, phi_min, phi_max); 


    if (i > one_percent + last_percent) 
    {
      printf("."); 
      fflush(stdout); 
      last_percent+=one_percent; 
    }

    if (event->should_save) t->Fill(); 
    if (event->pass) 
    {

      if (forced ) 
      {
        sum_weights += event->absorption_weight/event->lint; 
      }
      else
      {
        sum_weights += event->absorption_weight * event->interaction_weight; 

      }
    }

  }
  printf("\n"); 

  if (forced)
  {
    sum_weights *= sim.sampling_volume() / N ; 
  }
  else
  {
    sum_weights *= sim.sampling_area() / N;
  }

  t->Write(); 
  printf("Average Effective area over %gpi str is %g m^2\n", str, sum_weights); 
  TTree * ts = new TTree("meta","meta"); 
  const ToyDiscSetup * s = sim.setup(); 
  ts->Branch("setup", &s); 
  ts->Branch("forced", &forced); 
  ts->Branch("N", &N); 
  ts->Branch("mean_Aeff", &sum_weights); 
  ts->Fill(); 
  ts->Write(); 


  sim.getHSlant().Write(); 
}

