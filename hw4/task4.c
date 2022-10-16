double absoluteError(double actual, double measured) {
  return actual - measured;
}

double relativeError(double actual, double measured) {
  return absoluteError(actual, measured) / actual;
}
