class GMM {
  int K;
  int N;
  double[] pi;
  double[] mu;
  double[] va;
  GMM(int K, int N) {
    this.N = N;
    this.K = K;
    this.pi = new double[K];
    this.mu = new double[K];
    this.va = new double[K];
  }
  GMM(int K, int N, double[] pi, double[] mu, double [] va) {
    this.N = N;
    this.K = K;
    this.pi = new double[K];
    this.mu = new double[K];
    this.va = new double[K];
    this.setParams(pi, mu, va);
  }

  GMM(GMM g) {
    this.N = g.N;
    this.K = g.K;
    this.pi = new double[K];
    this.mu = new double[K];
    this.va = new double[K];
    for (int k = 0; k < K; k++) {
      this.pi[k] = g.pi[k];
      this.mu[k] = g.mu[k];
      this.va[k] = g.va[k];
    }
  }

  void setParams(double[] pi, double[] mu, double[] va) {
    this.pi = pi;
    this.mu = mu;
    this.va = va;
  }

  double sampleGMM() {
    int k = selectPi();
    double sd = Math.sqrt(this.va[k]);
    double sample = randomGaussian();
    return this.mu[k] + sample*sd;
  }

  double sampleGMM2() {
    int k = (int)random(0, K);
    double sd = Math.sqrt(this.va[k]);
    double sample = randomGaussian();
    return this.mu[k] + sample*sd;
  }

  double calculatePD(double x) {
    double ans = 0;
    for (int k = 0; k < K; k++) {
      double pd = (1/(Math.sqrt(2*Math.PI*this.va[k]))*Math.exp(-Math.pow(x-this.mu[k], 2)/(2*this.va[k])));
      ans += this.pi[k]*pd;
    }
    return ans;
  }

  int selectPi() {
    double[] cumPi = new double[K];
    cumPi[0] = pi[0];
    for (int i = 1; i < K; i++) {
      cumPi[i] = cumPi[i - 1] + pi[i];
    }

    float randomValue = random(1);

    for (int i = 0; i < K; i++) {
      if (randomValue < cumPi[i]) {
        return i;
      }
    }

    return K - 1;
  }

  void plotGMM(color c) {
    float dx = 0.1;
    stroke(c);
    strokeWeight(1);
    for (float x = 0; x < width; x += dx) {
      float X0 = map(x, 0, width, -3, 3);
      float Y0 = map((float)calculatePD(X0), 1, -0.5, 0, height);
      float X1 = map(x+dx, 0, width, -3, 3);
      float Y1 = map((float)calculatePD(X1), 1, -0.5, 0, height);
      line(x, Y0, x+dx, Y1);
    }
  }
}
