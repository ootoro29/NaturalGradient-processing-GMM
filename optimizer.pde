class Optimizer {
  GMM g;
  GMM preG;
  int K;
  int N;
  int step = 0;
  double X[];
  double[] param;
  double[] preParam;
  ArrayList<Float>P;
  boolean natural = true;
  boolean finish = false;
  Optimizer(GMM g, boolean natural) {
    this.g = new GMM(g);
    this.natural = natural;
    initOptimizer();
  }
  void initOptimizer() {
    preG = new GMM(g.K, g.N);

    //G.setPrams(dlist(0.6, 0.3, 0.2), dlist(-1, 0, 2), dlist(0.5, 0.1, 0.4));
    //preG.setPrams(dlist(0.8, 0.2, 0.1), dlist(-0.7, -0.2, 1), dlist(0.8, 0.2, 0.6));
    //preG.setPrams(dlist(0.1, 0.2, 0.7), dlist(-0.2, -0.5, 1.2), dlist(0.4, 0.4, 0.2));
    //preG.setPrams(dlist(1), dlist(-2.2), dlist(0.4));
    K = preG.K;
    N = preG.N;
    X = new double[N];
    param = new double[3*K];
    preParam = new double[3*K];
    for (int i = 0; i < K; i++) {
      param[i] = random(-2, 2);
    }
    for (int i = 0; i < K; i++) {
      param[K+i] = random(-2, 2);
    }
    for (int i = 0; i < K; i++) {
      param[2*K+i] = random(-2, 2);
    }
    setParam();
  }
  void setP(double[] p) {
    for (int i = 0; i < 3*K; i++) {
      param[i] = p[i];
    }
    setParam();
  }
  void setParam() {
    double[] mu = preG.mu;
    double[] va = preG.va;
    double[] pi = preG.pi;
    double sum = 0;
    for (int k = 0; k < K; k++) {
      mu[k] = param[k];
      va[k] = Math.exp(param[K+k]);
      pi[k] = Math.exp(param[2*K+k]);
      sum += pi[k];
    }
    for (int k = 0; k < K; k++) {
      pi[k] /= sum;
    }

    /*

     println("-----------------------");
     for (int k = 0; k < K; k++) {
     print(mu[k]+" ");
     }
     print("\n");
     for (int k = 0; k < K; k++) {
     print(va[k]+" ");
     }
     print("\n");
     for (int k = 0; k < K; k++) {
     print(pi[k]+" ");
     }
     print("\n");
     if (mu[0]!=mu[0])stop();
     println("-----------------------");
     */
  }

  Matrix DsoftMax() {
    Matrix softmax = new Matrix(K, 1);
    double sum = 0;
    for (int k = 0; k < K; k++) {
      softmax.Arr[k][0] = Math.exp(param[2*K+k]);
      sum += softmax.Arr[k][0];
    }
    for (int k = 0; k < K; k++) {
      softmax.Arr[k][0] /= sum;
    }
    Matrix ans = new Matrix(K, K);
    for (int i = 0; i < K; i++) {
      for (int j = 0; j < K; j++) {
        double A = 0;
        if (i == j) {
          A = softmax.Arr[i][0]*(1-softmax.Arr[i][0]);
        } else {
          A = -softmax.Arr[i][0]*softmax.Arr[j][0];
        }
        ans.Arr[i][j] = A;
      }
    }
    return ans;
  }
  void plotGMM(color c) {
    preG.plotGMM(c);
    noStroke();
    fill(255, 0, 0);
  }
  void update() {
    if (finish)return;
    P = new ArrayList();
    double[] mu = preG.mu;
    double[] va = preG.va;
    double[] pi = preG.pi;

    double X[] = new double[N];
    Matrix grad = new Matrix(3*K, 1);
    Matrix[] FN = new Matrix[N];
    Matrix[] GN = new Matrix[N];

    Matrix F = new Matrix(3*K, 3*K);
    if (step % 1 == 0) {
      for (int n = 0; n < N; n++) {
        X[n] = g.sampleGMM();
      }
    }

    for (int n = 0; n < N; n++) {
      FN[n] = new Matrix(3*K, 3*K);
      GN[n] = new Matrix(3*K, 1);
      Matrix gr = new Matrix(3*K, 1);
      double x = X[n];
      //P.add((float)x);
      double PD = preG.calculatePD(x);
      for (int k = 0; k < K; k++) {
        double exp = Math.exp(-(Math.pow(x-mu[k], 2))/(2*va[k]));
        //mu
        gr.Arr[k][0] = -pi[k]*(1/(Math.sqrt(2*Math.PI*va[k])))*((x-mu[k])/va[k])*exp;
        //va
        gr.Arr[K+k][0] = -pi[k]*(1/Math.sqrt(2*Math.PI))*exp*(((-1/(2*Math.pow(va[k], 1.5))))+((1/Math.sqrt(va[k]))*( Math.pow(x-mu[k], 2)/(2*Math.pow(va[k], 2)) )));
        //pi
        gr.Arr[2*K+k][0] = -(1/Math.sqrt(2*Math.PI*va[k]))*exp;

        gr.Arr[k][0] /= PD;
        gr.Arr[K+k][0] /= PD;
        gr.Arr[2*K+k][0] /= PD;
      }

      Matrix dLPi = new Matrix(K, 1);
      for (int k = 0; k < K; k++) {
        gr.Arr[K+k][0] *= Math.exp(param[K+k]);
        dLPi.Arr[k][0] = gr.Arr[2*K+k][0];
      }

      dLPi.equal_Mat(DsoftMax().mul_Mat(dLPi));

      for (int k = 0; k < K; k++) {
        gr.Arr[2*K+k][0] = dLPi.Arr[k][0];
      }

      FN[n].equal_Mat(gr.mul_Mat(gr.t_Mat()));
      GN[n].equal_Mat(gr);
    }

    for (int n = 0; n < N; n++) {
      for (int i = 0; i < 3*K; i++) {
        for (int j  =0; j < 3*K; j++) {
          F.Arr[i][j] += FN[n].Arr[i][j];
        }
        grad.Arr[i][0] += GN[n].Arr[i][0];
      }
    }

    for (int k = 0; k < K*3; k++) {
      for (int k2 = 0; k2 < K*3; k2++) {
        F.Arr[k][k2] /= N;
        if (k == k2)F.Arr[k][k2] += 0.001;
      }
    }

    for (int k = 0; k < K*3; k++) {
      grad.Arr[k][0] /= N;
    }
    Matrix IF = F.inv_Mat();


    if (natural)grad = IF.mul_Mat(grad);


    double alpha = (natural)?0.04:0.2;
    double alpha1 = (natural)?0.04:0.2;
    double alpha2 = (natural)?0.04:0.2;


    for (int k = 0; k < 3*K; k++) {
      preParam[k] = param[k];
    }

    for (int k = 0; k < K; k++) {
      param[k] -= alpha*grad.Arr[k][0];
      param[K+k] -= alpha1*grad.Arr[K+k][0];
      param[2*K+k] -= alpha2*grad.Arr[2*K+k][0];
    }

    Matrix par = new Matrix(3*K, 1);
    for (int k = 0; k < 3*K; k++) {
      par.Arr[k][0] = param[k]-preParam[k];
    }
    double sum = 0;
    sum = par.t_Mat().mul_Mat(F.mul_Mat(par)).Arr[0][0];
    double del = (natural)? 0.000006 : 0.000006;
    if (sum < del || step > 1000) {
      finish = true;
      if (natural) {
        println("自然勾配法:"+step+" "+dis(preG, g));
        saveHistoryG(+step+","+dis(preG, g));
      } else {
        println("確率勾配法:"+step+" "+dis(preG, g));
        saveHistoryN(+step+","+dis(preG, g));
      }
    }

    //println(sum, natural);

    step++;
    setParam();
  }

  double calP(double x) {
    double[] mu = preG.mu;
    double[] va = preG.va;
    double[] pi = preG.pi;
    double ans = 0;
    for (int k = 0; k < K; k++) {
      ans += pi[k]*(1/(Math.sqrt(2*Math.PI*va[k])))*(Math.exp(-Math.pow((x-mu[k]), 2)/(2*va[k])));
    }
    return ans;
  }

  double calQ(double x) {
    double[] mu = g.mu;
    double[] va = g.va;
    double[] pi = g.pi;
    double ans = 0;
    for (int k = 0; k < K; k++) {
      ans += pi[k]*(1/(Math.sqrt(2*Math.PI*va[k])))*(Math.exp(-Math.pow((x-mu[k]), 2)/(2*va[k])));
    }
    return ans;
  }
}
