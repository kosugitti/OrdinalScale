data{
  int L;      //Data Length
  int C;      //number of Categories
  int N;      //number of Subject
  int M;      //number of Items
  int D;      //number of dimensions
  int Pid[L]; //Index of Subject
  int Jid[L]; //Index of Items
  int<lower=1,upper=C> val[L];
}

parameters{
  vector<lower=0>[D] a[M];   //Discriminant Parameter
  ordered[C]         b[M];   //Difficulty Parameter
  vector[D]          theta[N];
}

transformed parameters{
  vector[C]  Pijk[N,M];
  for(n in 1:N){
    for(m in 1:M){
      for(c in 1:C){
        Pijk[n,m,c] = 0;
        for(c2 in 1:c){
          Pijk[n,m,c] = Pijk[n,m,c] + (a[m]' * theta[n] + b[m,c2]);
        }
      }
    }
  }

}

model{
  //likelihood
  for(l in 1:L){
    val[l] ~ categorical(softmax(Pijk[Pid[l],Jid[l]]));
  }
  //prior
  for(n in 1:N){
    for(d in 1:D){
      theta[n,d] ~ normal(0,1);
    }
  }
  for(m in 1:M){
    for(d in 1:D){
      a[m,d] ~ lognormal(0,sqrt(0.5));
    }
    for(c in 1:C){
      b[m,c] ~ normal(0,2);
    }
  }
}
