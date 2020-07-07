
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {

// data input

DATA_MATRIX(Lengths);
DATA_IMATRIX(surveyyear);
DATA_ARRAY(lambda);
DATA_VECTOR(Components);
DATA_VECTOR(agezero);

//parameters

PARAMETER(log_l); // mean of first observed cohort
PARAMETER(log_L); // mean of last observed cohort
PARAMETER(logit_k_reparam);//k for reparameterised vb for fournier & schnute
PARAMETER(log_sigma); // this and s1 are  used to give variable sd


//process parameters

Type l= exp(log_l);
Type L= exp(log_L) ;
Type k_reparam= invlogit(logit_k_reparam) ;
Type sigma= exp(log_sigma);


//std::cout << "lambda_dim" << std::endl << lambda.cols() << std::endl << std::endl;


//initialise vector of mus
int No_comp = Components.size() ;
int no_surveys = lambda.cols() ;

Type minagezero = min(agezero);

array <Type> mu_arr(1,No_comp,no_surveys);

//std::cout << "mu_dim" << std::endl << mu_arr.dim << std::endl << std::endl;


Type linf = (L-(l*(pow(k_reparam,No_comp-1))))/(1-(pow(k_reparam,No_comp-1)));
 //std::cout << "linf" << std::endl << linf << std::endl << std::endl;
  
Type K = -log(k_reparam);
   // std::cout << "K" << std::endl << K << std::endl << std::endl;
  
Type tzero =  minagezero-((1/log(k_reparam))*log((L-l)/(L-(l*(pow(k_reparam,No_comp-1))))));
    //std::cout << "tzero" << std::endl << tzero << std::endl << std::endl;


    for(int i =0;i<no_surveys;i++){
    
    mu_arr(0,0,i) = linf*(1-exp(-K*(agezero[i]-tzero)));
mu_arr(0,(No_comp-1),i) = linf*(1-exp(-K*((agezero[i]+(No_comp-1))-tzero)));
      
      
      for(int j=1; j<No_comp; j++){
        
        
        mu_arr(0,j,i) = mu_arr(0,j-1,i) + ((mu_arr(0,(No_comp-1),i) - mu_arr(0,0,i))*(((pow(k_reparam,(j+1)-2))-(pow(k_reparam,(j+1)-1)))/(1-(pow(k_reparam,No_comp-1)))));     //(j+1) for k_reparam^ as this effects calcs
      }  
      
      
    }
    
    
    
    
    
    
   // array <Type> mu_em1 = mu_arr.col(0);                  //uncomment to check mu_arr is being calculated properly
  //  std::cout << "mu_array" << std::endl << mu_em1.matrix() << std::endl << std::endl;
    
    
    
//no vector of sigmas needed




////////set up for nll fun

vector <Type> x =Lengths.col(2);
vector <Type> count =Lengths.col(3);

int n = x.size();

Type nll = 0.0; // initialize negative log likelihood
vector <Type> tmp(No_comp);

array <Type> lam1;
matrix <Type> lam2;
vector <Type> lam3;

int year = 0;
int survey =0;

for(int i = 0; i < n; i++){ 

year = surveyyear(i,1);
survey = surveyyear(i,0);

lam1 = lambda.col(survey);
lam2 = lam1.matrix();
lam3 = lam2.row(year); 


for(int j = 0; j < No_comp; j++){
  
  tmp[j] = ((lam3[j])*(dnorm(x[i], mu_arr(0,j,survey), sigma)));
  
}
nll -= count[i]*log(sum(tmp));
}
ADREPORT(l);
ADREPORT(L);
ADREPORT(k_reparam);
ADREPORT(sigma);
return nll;
}
