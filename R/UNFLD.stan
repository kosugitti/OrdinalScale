data{
  int L;      //Data Length
  int N;      //number of Subject
  int M;      //number of Items
  int D;      //number of dimensions
  int Pid[L]; //Index of Subject
  int Jid[L]; //Index of Items
  real val[L];
}

parameters{
  vector<lower=-5,upper=5>[D] a[M];   //coordinates
  vector[D] theta[N];
  real          beta0;
  real<lower=0> beta1;
  real<lower=0> sig;
}

transformed parameters{
  matrix[N,M] phi;
  for(n in 1:N){
    for(m in 1:M){
      phi[n,m] = dot_self(theta[n]-a[m])^0.5;
    }
  }
}

model{
  for(l in 1:L){
    val[l] ~ normal(beta0 - beta1*phi[Pid[l],Jid[l]] , sig);
  }

  for(n in 1:N){
    theta[n] ~ normal(0,1);
  }

  sig ~ student_t(4,0,5);
}
