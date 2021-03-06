
#include <TMB.hpp>

using namespace density;

template<class Type>
Type objective_function<Type>::operator() () {
  
  // data input
  
  DATA_MATRIX(Lengths);
  DATA_MATRIX(tau);
  DATA_IMATRIX(surveyyear);
  DATA_ARRAY(lambda);
  DATA_VECTOR(agezero);
  
  //parameters
  PARAMETER(log_l_mu); // 
  PARAMETER(log_L_mu); // 
  PARAMETER_MATRIX(lL);
  PARAMETER(raw_sdl);
  PARAMETER(raw_sdL);
  PARAMETER(raw_rho);  
  PARAMETER(logit_k_reparam);//k for reparameterised vb for fournier & schnute
  PARAMETER(log_s); // this and S are  used to give variable sd according to k_reparam
  PARAMETER(log_S);
  
  
  //process parameters
  
  vector <Type> l_par_vec = exp(log_l_mu + lL.col(0).array());
  vector <Type> L_par_vec = exp(log_L_mu + lL.col(1).array()) ;
  
  Type sdl = exp(raw_sdl);
  Type sdL = exp(raw_sdL);
  Type rho = Type(-1.0) + Type(2.0) * invlogit(raw_rho); // (-1,1)
  Type k_reparam= invlogit(logit_k_reparam) ;
  Type s= exp(log_s);
  Type S= exp(log_S);
  
 //std::cout << "l par vec" << std::endl << l_par_vec << std::endl << std::endl;
  //std::cout << "L par vec" << std::endl << L_par_vec << std::endl << std::endl;
  
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
  
  Type K = -log(k_reparam);
  // std::cout << "K" << std::endl << K << std::endl << std::endl;
  
  Type tzero ;
  //std::cout << "tzero" << std::endl << tzero << std::endl << std::endl;
  
  
  for(int i =0; i < (No_comp-1);i++){
    mu(i,0,0) = l_par_vec.sum()/l_par_vec.size();
    mu(i+(No_comp-1),No_comp-1,0) = L_par_vec.sum()/L_par_vec.size();
  }
  
  
  for(int i = 0; i < No_years;i++){
    mu(i+(No_comp-1),0,0) = l_par_vec[i] ;
    mu(i+(2*(No_comp-1)),No_comp-1,0) = L_par_vec[i];
  }
  
  if(no_surveys >1){
    for(int j = 1; j<no_surveys;j++){
      for(int i = 0; i < No_maxyearbff;i++){
        linf = (mu(i+No_comp-1,No_comp-1,0)-(mu(i,0,0)*(pow(k_reparam,No_comp-1))))/(1-(pow(k_reparam,No_comp-1)));
        tzero =  agezero[0]-((1/log(k_reparam))*log((mu(i+No_comp-1,No_comp-1,0)-mu(i,0,0))/(mu(i+No_comp-1,No_comp-1,0)-(mu(i,0,0)*(pow(k_reparam,No_comp-1))))));
        
        mu(i,0,j) = linf*(1-exp(-K*(agezero[j]-tzero)));
        mu(i+No_comp-1,No_comp-1,j) =linf*(1-exp(-K*((agezero[j]+(No_comp-1))-tzero))); 
      }
    } 
  }
    
    //mu.block((No_comp-1),0,No_years,1) = l_par_vec;		
    
    //mu.block(No_comp-1,No_comp-1,No_years,1) = L_par_vec;	//assign last column from No_comp-1 to No_mycomp with L vector
    
    //set up diffiw, l and L for strorage in loop
    
    int diffiw=0;
  Type l=0;
  Type L=0;
  for(int j =0; j <no_surveys; j++){
    for(int i=1;i<No_mycomp;i++){
      
      for(int w=1;w<(No_comp-1);w++){
        
        diffiw = i-w;
        if(diffiw >= 0 && diffiw < No_maxyearbff){
          l=mu(diffiw,0,j);
          L=mu(diffiw + No_comp - 1, No_comp -1,j);
        }
        else {
          l=0;		//zero so should flag up if code not working properly
          L=0;
        }
        
        mu(i,w,j)= mu(i-1,w-1,j) + ((L - l)*(((pow(k_reparam,(w+1)-2))-(pow(k_reparam,(w+1)-1)))/(1-(pow(k_reparam,No_comp-1)))));
        
      }
    }
  }
  
  //std::cout << "mu" << std::endl << mu.col(2).matrix() << std::endl << std::endl;
  
  int startrow = No_comp-1;
  
  //matrix <Type> new_mu = mu.block(startrow, 0, No_years, No_comp);	
  
  //std::cout << "new mu" << std::endl << new_mu << std::endl << std::endl;
  
  
  //calculate array of sigmas as well
  array <Type> SD(No_mycomp,No_comp,no_surveys);
  
  for(int j =0; j <no_surveys; j++){
    for(int i=No_comp-1;i<No_maxyearbff;i++){
      for(int w=0;w<(No_comp);w++){
        
        SD(i,w,j) = s+ (S-s)*((mu(i,w,j) - exp(log_l_mu))/(exp(log_L_mu) - exp(log_l_mu))) ;
      }
    }
  }
  
  //std::cout << "sd" << std::endl << SD.col(1).matrix() << std::endl << std::endl;
  
  
  ////////set up for nll fun
  
  
  
  vector <Type> x =Lengths.col(2);
  vector <Type> count =Lengths.col(3);
  
  int n = x.size(); 
  int t = No_comp;
  
  Type nll = 0.0; // initialize negative log likelihood
  
  // random effects component
  matrix<Type> Sigma(2,2);
  Sigma(0,0) = pow(sdl, 2.0);
  Sigma(1,1) = pow(sdL, 2.0);
  Sigma(1,0) = rho * sdl * sdL;
  Sigma(0,1) = rho * sdl * sdL;
  for(int j = 0; j < No_years; j++){
    nll += MVNORM(Sigma)(vector<Type>(lL.row(j))); // note + as return negative log density
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
      
      
      nll -= count[i]*(tau(i,j)*((log(lam3[j])) + (dnorm(x[i], new_mu2(year,j), new_SD2(year,j), true))));
    }
  }
  ADREPORT(log_l_mu + lL.col(0).array());
  ADREPORT(log_L_mu + lL.col(1).array());
  return nll;
}

