#include <iostream>
#include <math.h>
using namespace std;

///////////////////////////////////////////////////////////////////////////
//
// Class: Polynomial
//

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

class Polynomial {
public:
  Polynomial dx() const;                                // Derivative of Polynomial
  
  float operator()(const float x) const;                // Evaluate Polynomial
  Polynomial operator+(const Polynomial& p) const;      // Add 
  Polynomial operator*(const Polynomial& p) const;      // Multiply
  void operator=(const Polynomial& p);                  // Assignment
  bool operator==(const Polynomial& p) const;           // Equality (relative, not exact, since float coeff)
  
  Polynomial(const Polynomial& p);                      // Copy constructor
  Polynomial(const int degree, const float coeff[]);    // Create from an array (coeff from An to A0)
  ~Polynomial();                                        // Destructor: get rid of the coeff array

  friend std::ostream& operator<<(std::ostream& os, const Polynomial& p);  // Polynomial printing

private:
  // Come from struct and hacked quickly; we don't want to support them, though
  // If you change anything in "private" you will need to submit ClassyPoly.h
  // as well as ClassyPoly.cpp by combining them in a zip file and submitting
  // that
  // Also, if you change the way the data is stored, you will likely want to
  // modify the way the operator<<() function works
  //
  float evaluate(const float x) const;                  // Evaluate Polynomial
  Polynomial add(const Polynomial& p) const;            // Add (from struct)
  Polynomial derivative() const;                        // Derivative (from struct) 

  
  int    _degree;
  float* _coeff;
};


 std::ostream& operator<<(std::ostream& os, const Polynomial& p) {
  if (p._degree < 0) {
    std::cerr << "Error: Attempted to output deleted polynomial; exiting";
    exit(-1);
  }
  for (int i = 0; i < p._degree; ++i) {
    if (p._coeff[i] == 1) 
      os << "x";
    else
      os << p._coeff[i] << "x";
    if ((p._degree - i) > 1) 
      os <<  "^" <<  (p._degree - i);
    os << " + ";
  }
  os << p._coeff[p._degree];
  return os;
}

#endif


///////////////////////////////////////////////////////////////////////////
//
// Public operators and methods

Polynomial Polynomial::dx() const { //derivative of polynomial
  return Polynomial::derivative();
}

float Polynomial::operator()(const float x) const { //evaluate polynomial
  return Polynomial::evaluate(x);
}

Polynomial Polynomial::operator + (const Polynomial & p) const { //add polynomial
  return Polynomial::add(p);
}

Polynomial Polynomial::operator * (const Polynomial & p) const { //multiply polynomial
  if (_degree == 0) { //if the degree is 0, the second polynomial is a constant, so we have to add it to itself
    Polynomial product(p);

    for (int i = 1; i < _coeff[0]; i++) {
      product = product.add(p);
    }

    return product;
  }

  if (p._degree == 0) { //if the degree is 0, the second polynomial is a constant, so we have to add it to itself
    Polynomial product( * this);

    for (int i = 1; i < p._coeff[0]; i++) {
      product = product.add( * this);
    }

    return product;
  }

  int productDegree = p._degree + _degree;
  float * productCoeff = new float[productDegree + 1];

  for (int i = 0; i <= productDegree; i++) { // fill it with zeros
    productCoeff[i] = 0;
  }

  if (_degree > p._degree) {
    for (int i = 0; i <= _degree; i++) {
      for (int j = 0; j <= p._degree; j++) {
        productCoeff[i + j] += _coeff[i] * p._coeff[j];
      }
    }
  } else {
    for (int i = 0; i <= p._degree; i++) {
      for (int j = 0; j <= _degree; j++) {
        productCoeff[i + j] += _coeff[j] * p._coeff[i];
      }
    }
  }

  Polynomial product(productDegree, productCoeff);
  return product;
}

void Polynomial::operator = (const Polynomial & p) { //assignment
  for (int i = 0; i <= _degree; i++) {
    _coeff[i] = 0;
  }

  _degree = 0;
  delete[] _coeff;

  Polynomial pNew(p);

  _degree = p._degree;
  _coeff = new float[_degree + 1];

  for (int i = 0; i <= pNew._degree; i++) {
    _coeff[i] = pNew._coeff[i];
  }
}

bool Polynomial::operator == (const Polynomial & p) const { //equality
  if (_degree < 0 || p._degree < 0) {
    return false;
  }

  if (_degree != p._degree) { //if the degrees of both polynomials are different, they can't be equal
    return false;
  }

  if (_degree == p._degree) {
    for (int i = 0; i <= _degree; i++) { //loop thru every coefficient til the term of the highest degree
      if (fabs(_coeff[i] - p._coeff[i]) > 0.00001) { //if the coefficients aren't the same, the 2 polynomials aren't equal
        return false;
      }
    }
  }

  return true; //if you reach this point, then the polynomials must have all the same coefficients and same degree
}

///////////////////////////////////////////////////////////////////////////
//
// Constructors and Destructor
//
// Note: no default constructor
//

Polynomial::Polynomial(const Polynomial & p) { //copy constructor
  _degree = p._degree;
  _coeff = new float[p._degree + 1];

  for (int i = 0; i <= p._degree; i++) {
    _coeff[i] = p._coeff[i];
  }

}

Polynomial::Polynomial(const int degree,
  const float coeff[]) { //create from an array (coeff from An to A0)
  if (degree < 0) { //if this polynomial is invalid, make a new even more invalid polynomial
    _degree = degree;
    _coeff = new float[1];
    _coeff[0] = 0;
    return;
  }

  _degree = degree;
  _coeff = new float[_degree + 1];

  for (int i = 0; i <= degree; i++) {
    _coeff[i] = coeff[i];
  }
}

Polynomial::~Polynomial() { //destructor
  for (int i = 0; i <= _degree; i++) {
    _coeff[i] = 0;
  }

  _degree = 0;
  delete[] _coeff;
}

///////////////////////////////////////////////////////////////////////////
//
// Private Methods
//

float Polynomial::evaluate(const float x) const {
  float total = 0;

  if (_degree < 0) {
    return NAN;
  }

  int i = 0;
  while (i <= _degree) {
    total = x * (total) + _coeff[i];
    i++;
  }

  return total;
}

Polynomial Polynomial::add(const Polynomial & p) const {
  if (_degree < 0 || p._degree < 0) {
    float * invalidCoeff = new float[1];
    invalidCoeff[0] = 0;
    Polynomial invalid(-1, invalidCoeff);
    return invalid;
  }

  Polynomial sum( * this); //if _degree >= p._degree
  Polynomial sumAlt(p); //if _degree < p._degree

  int i = 0;
  int diff = 0;
  int P1Bigger = 1;

  if (_degree > p._degree) {
    diff = _degree - p._degree;
  }

  if (_degree < p._degree) {
    diff = p._degree - _degree;
    P1Bigger = 0;
  }

  if (diff > 0) {
    while (i < diff) {
      if (P1Bigger) {
        sum._coeff[i] = p._coeff[i];
      }

      if (!P1Bigger) {
        sumAlt._coeff[i] = _coeff[i];
      }
      i++;
    }
  }

  if (P1Bigger) {
    int j = 0;
    for (int i = diff; i <= _degree; i++) {
      sum._coeff[i] = _coeff[i] + p._coeff[j];
      j++;
    }
    return sum;
  }

  if (!P1Bigger) {
    int j = 0;
    for (int i = diff; i <= p._degree; i++) {
      sumAlt._coeff[i] = _coeff[j] + p._coeff[i];
      j++;
    }
    return sumAlt;
  }

  return sum;
}

Polynomial Polynomial::derivative() const {
  /*if (_degree < 0 || p._degree < 0)  {
    _coeff = new float [0];
    _degree = -1;
    return Polynomial;
  }

  if (_degree == 0) {
    _coeff = new float [1];
    _coeff[0] = 0;
    _degree = 0;
    return Polynomial;
  }
  */

  Polynomial derive( * this);

  derive._degree = _degree - 1;
  derive._coeff = new float[_degree + 1];

  int i = 0;
  while (i <= _degree) {
    derive._coeff[i] = _coeff[i] * (_degree - i);
    i++;
  }

  return derive;
}

///////////////////////////////////////////////////////////////////////////
//
// Test driver
// Some very limited testing; should test ==
// 

#ifndef MARMOSET_TESTING

int main() {
  float coeff[] = {
    1,
    2,
    3,
    4
  }; // x^2 + 2x + 3
  Polynomial p0(1, coeff);
  Polynomial p1(2, coeff);
  Polynomial p2(3, coeff);

  cout << "When x = 2.2, \"" << p0 << "\" evalates to: " << p0(2.2) << endl << endl;
  cout << "When x = 2.2, \"" << p1 << "\" evalates to: " << p1(2.2) << endl << endl;
  cout << "When x = 2.2, \"" << p2 << "\" evalates to: " << p2(2.2) << endl << endl;

  Polynomial p = p1 + p2;
  Polynomial q = p1 * p2;
  cout << "p1 + p2 =  " << p << endl << endl;
  cout << "p1 * p2 =  " << q << endl << endl;

  cout << "dp/dx =  " << p.dx() << endl;
  cout << "dq/dx =  " << q.dx() << endl;

  return 0;
}

#endif