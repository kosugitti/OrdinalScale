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
  vector[D] raw_a[M-3];   //coordinates
  real<lower=0> a11;
  real<lower=0> a12;
  real<upper=0> a21;
  real<upper=0> a22;
  real<lower=0> a31;
  real<upper=0> a32;
  vector[D] theta[N];
  real          beta0;
  real<lower=0> beta1;
  real<lower=0> sig;
}

transformed parameters{
  matrix[N,M] phi;
  vector[D] a[M];

  //dim constraint
  for(m in 1:(M-3)){
    a[m] = raw_a[m];
  }
  a[M,1] = a11;
  a[M,2] = a12;
  a[(M-1),1] = a21;
  a[(M-1),2] = a22;
  a[(M-2),1] = a31;
  a[(M-2),2] = a32;

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

  for(m in 1:(M-3)){
    raw_a[m] ~ normal(0,1);
  }
  a11 ~ normal(0,1);
  a12 ~ normal(0,1);
  a21 ~ normal(0,1);
  a22 ~ normal(0,1);
  a31 ~ normal(0,1);
  a32 ~ normal(0,1);

  for(n in 1:N){
    theta[n] ~ normal(0,1);
  }
  sig ~ student_t(4,0,5);
}
