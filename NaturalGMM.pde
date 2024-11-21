GMM G;
Optimizer O;
Optimizer O2;
void setup() {
  size(600, 600);
  //G.setPrams(dlist(1), dlist(-1), dlist(0.3));
  //G = new GMM(3, 1000, dlist(0.6, 0.3, 0.1), dlist(-2, 0, 1), dlist(0.2, 0.9, 0.6));
  G = new GMM(4, 1000, dlist(0.4, 0.1, 0.2, 0.3), dlist(-2, 1.5, 0, 1), dlist(0.2, 1.2, 0.9, 0.1));

  //G = new GMM(8, 1000, dlist(0.4, 0.05, 0.12, 0.03, 0.1, 0.08, 0.12, 0.1), dlist(-2, 1.5, 0, 1, 3, -0.25, -0.6, -1.8), dlist(0.2, 1.2, 0.9, 0.1, 0.4, 0.2, 0.05, 3));

  //G.setPrams(dlist(0.6, 0.3, 0.2), dlist(-1, 0, 2), dlist(0.1, 0.1, 0.4));
  //G = new GMM(4, 100000);
  //G.setPrams(dlist(0.4, 0.15, 0.25, 0.2), dlist(-1, 0.5, 0, 2), dlist(0.1, 1.2, 0.1, 0.4));
  //G.setPrams(dlist(1), dlist(-1), dlist(0.5));
  optInit();
}
void optInit() {
  O = new Optimizer(G, true);
  O2 = new Optimizer(G, false);

  double P[] = new double[30];
  for (int i = 0; i < 30; i++) {
    P[i] = random(-2, 2);
  }
  O.setP(P);
  O2.setP(P);
}
void draw() {
  background(255);
  G.plotGMM(color(255, 0, 0));
  O.plotGMM(color(0, 0, 255));
  O2.plotGMM(color(0, 255, 0));
  for (int i = 0; i < 1; i++) {
    O.update();
    O2.update();
  }
  if (O.finish && O2.finish) {
    optInit();
  }
}
