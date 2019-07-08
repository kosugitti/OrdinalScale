data{
  int P;  // Number of items
  int K;  // Number of response Categories
  int N;  // Number of subjects
  int<lower=1,upper=K> y[N,P]; // responses
}

parameters{
  real<lower=0,upper=5> a[P]; // discriminant
  ordered[K-1] b[P];          // difficulty
  vector<lower=-5,upper=5>[N] theta; // latant trait
}

transformed parameters{
  real b[6,5];
  vector<lower=0,upper=1>[5-1] pa[Ni,6];
  simplex[5] p[Ni,6];

  for(j in 1:6){
    for(c in 1:5){
      if(c==1){
        b[j,c] = ba[j,c];
      }else if(c == 5){
        b[j,c] = ba[j,c-1];
      }else{
        b[j,c] = (ba[j,c-1]+ba[j,c])/2;
      }
    }
  }

  for(i in 1:Ni){
    for(j in 1:6){
      for(c in 1:4){
        pa[i,j,c] = 1/(1+exp(-1.7*a[j]*(theta[i]-ba[j,c])));
      }
    }
  }

  for(i in 1:Ni){
    for(j in 1:6){
      for(c in 1:5){
        if(c==1){
          p[i,j,c] = 1-pa[i,j,c];
        }else if(c==5){
          p[i,j,c] =  pa[i,j,c-1];
        }else{
          p[i,j,c] = pa[i,j,c-1] - pa[i,j,c];
        }
      }
    }
  }

}

model{
  for(i in 1:Ni){
    theta[i] ~ normal(0,1);
  }

  for (j in 1:6){
		a[j]~lognormal(0,sqrt(0.5));
		for (c in 1:4){
			ba[j,c]~normal(0,2);
		}
	}

	for(l in 1:L){
	  Val[l] ~ categorical(p[Pid[L],Iid[l]]);
	}
}


generated quantities{
  real log_lik[L];
  real pred[L];
  for(l in 1:L){
	  log_lik[l] = categorical_lpmf(Val[l]|p[Pid[L],Iid[l]]);
	  pred[l] = categorical_rng(p[Pid[L],Iid[l]]);
	}

}
