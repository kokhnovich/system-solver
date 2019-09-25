//
// Created by user on 3.09.19.
//
#pragma once

#ifndef SYSTEM_SOLVER__FRACTION_H_
#define SYSTEM_SOLVER__FRACTION_H_

long long lcd(long long a, long long b) {
  return a / __gcd(a, b) * b;
}

class Fraction {
 public:
  long long numerator, denominator;

  Fraction(long long n = 1, long long d = 1) : numerator(n), denominator(d) {
    if (denominator == 0) throw logic_error("divide by zero");
    Normalize();
  }

  void Normalize() {
    if (1ll * numerator * denominator < 0) {
      numerator = -1 * abs(numerator);
      denominator = abs(denominator);
    } else {
      numerator = abs(numerator);
      denominator = abs(denominator);
    }
    long long gcd = __gcd(abs(numerator), abs(denominator));
    numerator /= gcd;
    denominator /= gcd;
  }

  Fraction operator=(const Fraction& other) {
    numerator = other.numerator, denominator = other.denominator;
    return {numerator, denominator};
  }

  Fraction operator+(const Fraction& otherFraction) const {
    long long n = numerator * otherFraction.denominator + otherFraction.numerator * denominator;
    long long d = denominator * otherFraction.denominator;
    return {n, d};
  }
  Fraction operator-(const Fraction& otherFraction) const {
    long long n = numerator * otherFraction.denominator - otherFraction.numerator * denominator;
    long long d = denominator * otherFraction.denominator;

    return {n, d};
  }
  Fraction operator*(const Fraction& otherFraction) const {
    long long n = numerator * otherFraction.numerator;
    long long d = denominator * otherFraction.denominator;
    return {n, d};
  }
  Fraction operator/(const Fraction& otherFraction) const {
    long long n = numerator * otherFraction.denominator;
    long long d = denominator * otherFraction.numerator;
    return {n, d};
  }

  bool operator==(const Fraction& otherFraction) const {
    return numerator == otherFraction.numerator && denominator == otherFraction.denominator;
  }
  bool operator!=(const Fraction& otherFraction) const {
    return !(*this == otherFraction);
  }
  bool operator<(const Fraction& other) const {
    long long lcd_ = lcd(denominator, other.denominator);
    return numerator * (lcd_ / denominator) < other.numerator * (lcd_ / denominator);
  }
  bool operator>(const Fraction& other) const {
    return !(*this < other) && !(*this == other);
  }

  Fraction operator+=(const Fraction& otherFraction) {
    *this = *this + otherFraction;
    return *this;
  }
  Fraction operator-=(const Fraction& otherFraction) {
    *this = *this - otherFraction;
    return *this;
  }
  Fraction operator/=(const Fraction& otherFraction) {
    *this = *this / otherFraction;
    return *this;
  }

  ostream& operator<<(ostream& os) const {
    os << setw(7) << numerator << "/" << denominator << endl;
    return os;
  }

  operator long long() const {
    return *this;
  }
  operator int() const {
    return *this;
  }

  void show() {
    cout << setw(7) << numerator << "/" << denominator << endl;
  }
};

ostream& operator<<(ostream& os, const Fraction& dt) {
  // os << fixed << setprecision(5) << (double) dt.numerator / dt.denominator;
  os << setw(5) << dt.numerator << '/' << dt.denominator;
  return os;
}

Fraction abs(const Fraction& val) {
  return Fraction(abs(val.numerator), abs(val.denominator));
}


#endif //SYSTEM_SOLVER__FRACTION_H_
