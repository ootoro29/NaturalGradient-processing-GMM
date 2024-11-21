double[] dlist(double ...a) {
  return a;
}

double dis(GMM G1, GMM G2) {
  double del = 0.0001;
  double ans = 0;
  for (double x = -10; x < 10; x +=del) {
    ans += del*(G1.calculatePD(x)*Math.log(G1.calculatePD(x)/G2.calculatePD(x)));
  }
  return ans;
}

String folderName = "res/4";
ArrayList<String>historyN = new ArrayList();

void saveHistoryN(String s) {
  historyN.add(s);
  String[] Save = new String[historyN.size()];
  for (int i = 0; i < historyN.size(); i++) {
    Save[i] = historyN.get(i);
  }
  saveStrings(folderName+"/historyN.txt", Save);
}

ArrayList<String>historyG = new ArrayList();

void saveHistoryG(String s) {
  historyG.add(s);
  String[] Save = new String[historyG.size()];
  for (int i = 0; i < historyG.size(); i++) {
    Save[i] = historyG.get(i);
  }
  saveStrings(folderName+"/historyG.txt", Save);
}
