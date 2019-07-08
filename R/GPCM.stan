data{
  int P;  // Number of items
  int K;  // Number of response Categories
  int N;  // Number of subjects
  int<lower=1,upper=K> y[N,P]; // responses 
}

parameters{
  real theta[N];    // latant trait
  vector[K-1] b[P];  // difficulty
  real<lower=0> a[P]; //discriminant
}

transformed parameters{
  vector[K] b2[P];
  for(p in 1:P){
    b2[p,1] = 0;
    for(k in 2:K){
      b2[p,k] = b[p,(k-1)];
    }
  }
}

model{
  vector[K] pi;
  for(n in 1:N){
    for(p in 1:P){
      pi = rep_vector(0,K);
      for(k in 1:K){
        for(k2 in 1:k){
          pi[k] = pi[k] + (a[p] * (theta[n] - b2[p,k2]));
        }
      }

      y[n,p] ~ categorical(softmax(pi));
    }
  }

  // prior
  theta ~ normal(0,1);
  a ~ student_t(4,0,10);
  for(k in 1:(K-1)){
    for(p in 1:P){
      b[p,k] ~ normal(0,100);
    }
  }
}

