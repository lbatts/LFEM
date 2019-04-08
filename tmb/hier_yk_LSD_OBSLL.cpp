#include <TMB.hpp>

using namespace density;

template<class Type>
Type objective_function<Type>::operator() () {
  
  // data input
  
  DATA_MATRIX(Lengths);
  DATA_IMATRIX(surveyyear);
  DATA_ARRAY(lambda);
  DATA_VECTOR(agezero);
  
  //parameters
  PARAMETER(log_l_mu); // 
  PARAMETER(log_L); // 
  PARAMETER_MATRIX(lk);
  PARAMETER(raw_sdl);
  PARAMETER(raw_sdk);
  PARAMETER(raw_rho);  
  PARAMETER(logit_k_reparam_mu);//k for reparameterised vb for fournier & schnute
  PARAMETER(log_s); // this and S are  used to give variable sd according to k_reparam
  PARAMETER(log_S);
  
  
  
  //process parameters
  
  vector <Type> l_par_vec = exp(log_l_mu + lk.col(0).array());
  Type L = exp(log_L);
  Type l_mu = exp(log_l_mu);
  vector <Type> k_par_vec_2 = logit_k_reparam_mu + lk.col(1).array();	//needs to be initialised this way for some reason
  vector <Type> k_par_vec = invlogit(k_par_vec_2) ;
  Type sdl = exp(raw_sdl);
  Type sdk = exp(raw_sdk);
  Type rho = Type(-1.0) + Type(2.0) * invlogit(raw_rho); // (-1,1)
  Type s= exp(log_s);
  Type S= exp(log_S);
  
  
  //std::cout << "l par vec" << std::endl << l_par_vec << std::endl << std::endl;
  //std::cout << "k par vec" << std::endl << k_par_vec << std::endl << std::endl;
  
  //initialise integers for manipulation etc.
  int No_comp = lambda.dim(1);
  int No_years = lambda.dim(0);
  int no_surveys = lambda.cols() ;
  int No_maxyearbff = No_years + No_comp - 1;
  int No_mycomp = No_maxyearbff + No_comp -1;
  
  
  //matrix <Type> mu(No_mycomp,No_comp);	//set up mu matrix which is all years + 2*(Number of components -1)
  array <Type> mu(No_mycomp,No_comp,no_surveys);
  
  //std::cout << "mu_dim" << std::endl << mu.dim << std::endl << std::endl;
  
  Type linf; 
  //std::cout << "linf" << std::endl << linf << std::endl << std::endl;
  
  Type K;//     
  // std::cout << "K" << std::endl << K << std::endl << std::endl;
  
  Type tzero ;
  //std::cout << "tzero" << std::endl << tzero << std::endl << std::endl;
  
  
  //vector <Type> k_rep_full(No_mycomp-1);
  //k_rep_full.head(No_comp-1) = k_par_vec.sum()/k_par_vec.size();//mean of random effects
  //k_rep_full.tail(No_comp-1) = k_par_vec.sum()/k_par_vec.size();
  //k_rep_full.segment(No_comp-1,No_years) = k_par_vec;
  //std::cout << "k rep full" << std::endl << k_rep_full << std::endl << std::endl;
  
  
  
  vector <Type> k_par_vec_plus(No_mycomp);
  
  k_par_vec_plus.tail(No_comp-1) = k_par_vec.sum()/k_par_vec.size();//mean of random effects
  
  for(int i =0; i < (No_comp-1);i++){
    mu(i,0,0) = l_par_vec.sum()/l_par_vec.size();
    mu(i+(No_comp-1),No_comp-1,0) = L;
    k_par_vec_plus[i]=k_par_vec.sum()/k_par_vec.size();
  }
  
  
  for(int i = 0; i < No_years;i++){
    mu(i+(No_comp-1),0,0) = l_par_vec[i] ;
    mu(i+(2*(No_comp-1)),No_comp-1,0) = L;
    k_par_vec_plus[i+(No_comp-1)] = k_par_vec[i];
  }
  
  //std::cout << "k par vec plus" << std::endl << k_par_vec_plus << std::endl << std::endl;
  
  for(int j = 1; j<no_surveys;j++)
    for(int i = 0; i < No_maxyearbff;i++){
      K = -log((k_par_vec_plus.segment(i,No_comp).sum())/No_comp);
      linf = (mu(i+No_comp-1,No_comp-1,0)-(mu(i,0,0)*(pow((k_par_vec_plus.segment(i,No_comp).sum())/No_comp,No_comp-1))))/(1-(pow((k_par_vec_plus.segment(i,No_comp).sum())/No_comp,No_comp-1)));
      tzero =  agezero[0]-((1/log((k_par_vec_plus.segment(i,No_comp).sum())/No_comp))*log((mu(i+No_comp-1,No_comp-1,0)-mu(i,0,0))/(mu(i+No_comp-1,No_comp-1,0)-(mu(i,0,0)*(pow((k_par_vec_plus.segment(i,No_comp).sum())/No_comp,No_comp-1))))));
      
      //std::cout << "cohort k rep mean" << std::endl << ((k_par_vec_plus.segment(i,No_comp).sum())/No_comp) << std::endl << std::endl;
      
      mu(i,0,j) = linf*(1-exp(-K*(agezero[j]-tzero)));
      mu(i+No_comp-1,No_comp-1,j) =linf*(1-exp(-K*((agezero[j]+(No_comp-1))-tzero))); 
    }
    
    //std::cout << "mu" << std::endl << mu.col(1).matrix() << std::endl << std::endl;
  //std::cout << "k par vec plus" << std::endl << k_par_vec_plus << std::endl << std::endl;
  
  
  
  //set up diffiw, l and k_reparam for strorage in loop
  int diffiw=0;
  Type l=0;
  Type L_tmp=0;
  Type k_reparam = 0;
  
  for(int j =0; j <no_surveys; j++){
    for(int i=1;i<No_mycomp;i++){
      
      for(int w=1;w<(No_comp);w++){
        
        diffiw = i-w;
        if(diffiw >= 0 && diffiw < No_maxyearbff){
          l = mu(diffiw,0,j);
          L_tmp = mu(diffiw + No_comp - 1, No_comp -1,j);
          k_reparam = k_par_vec_plus[i-1];
        }
        else {
          l=0;		//zero so should flag up if code not working properly
          L_tmp=0;
          k_reparam = 0;
        }
        
        mu(i,w,j)= mu(i-1,w-1,j) + ((L_tmp - l)*(((pow(k_reparam,(w+1)-2))-(pow(k_reparam,(w+1)-1)))/(1-(pow(k_reparam,No_comp-1)))));
        
      }
    }
  }
  
  //std::cout << "mu" << std::endl << mu.col(2).matrix() << std::endl << std::endl;
  
  int startrow = No_comp-1;
  
  
  //calculate array of sigmas as well
  array <Type> SD(No_mycomp,No_comp,no_surveys);
  
  for(int j =0; j <no_surveys; j++){
    for(int i=No_comp-1;i<No_maxyearbff;i++){
      for(int w=0;w<(No_comp);w++){
        
        SD(i,w,j) = s+ (S-s)*((mu(i,w,j) - exp(log_l_mu))/(L - exp(log_l_mu))) ;
      }
    }
  }
  //std::cout << "sd" << std::endl << SD.col(2).matrix() << std::endl << std::endl;
  
  ////////set up for nll fun
  
  
  
  vector <Type> x =Lengths.col(2);
  vector <Type> count =Lengths.col(3);
  
  int n = x.size(); 
  int t = No_comp;
  
  Type nll = 0.0; // initialize negative log likelihood
  vector <Type> tmp(t);
  
  // random effects component
  matrix<Type> Sigma(2,2);
  Sigma(0,0) = pow(sdl, 2.0);
  Sigma(1,1) = pow(sdk, 2.0);
  Sigma(1,0) = rho * sdl * sdk;
  Sigma(0,1) = rho * sdl * sdk;
  for(int j = 0; j < No_years; j++){
    nll += MVNORM(Sigma)(vector<Type>(lk.row(j))); // note + as return negative log density
  }
  
  
  
  ////length class likelihood component
  
  array <Type> lam1;
  matrix <Type> lam2;
  vector <Type> lam3;
  array <Type> new_mu;
  matrix <Type> new_mu2;
  array <Type> new_SD;
  matrix <Type> new_SD2;
  
  
  int year = 0;
  int survey =0;
  
  for(int i = 0; i < n; i++){ 
    
    year = surveyyear(i,1);
    survey = surveyyear(i,0);
    
    lam1 = lambda.col(survey);
    lam2 = lam1.matrix();
    lam3 = lam2.row(year); 
    new_mu = mu.col(survey);
    new_mu2 = new_mu.matrix().block(startrow, 0, No_years, No_comp);
    new_SD = SD.col(survey);
    new_SD2 = new_SD.matrix().block(startrow, 0, No_years, No_comp);
    
    for(int j = 0; j < t; j++){
      
      tmp[j] = ((lam3[j])*(dnorm(x[i], new_mu2(year,j), new_SD2(year,j))));
      
    }
    nll -= count[i]*log(sum(tmp));	//is a plus after the nll and not a minus so returns negative log lik
  }
  //std::cout << "lambda" << lambda << std::endl << std::endl;
  //std::cout << "x" << x << std::endl << std::endl;
  //std::cout << "new_mu" << new_mu << std::endl << std::endl;
  //std::cout << "comp_sigma" << comp_sigma << std::endl << std::endl;
  //std::cout << "count" << count << std::endl << std::endl;
  ADREPORT(k_par_vec_2);
  ADREPORT(log_l_mu + lk.col(0).array());
  ADREPORT(l_mu);
  ADREPORT(L);
  ADREPORT(sdl);
  ADREPORT(sdk); 
  ADREPORT(k_reparam);
  ADREPORT(S);
  ADREPORT(s);
  return nll;
}

