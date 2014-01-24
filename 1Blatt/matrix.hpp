#ifndef MATRIX_HPP
#define MATRIX_HPP

class Matrix {
private:
  const unsigned n_;
public:
  // constructor
  Matrix() : n_(0) { }
  Matrix(const unsigned n) : n_(n) { }

  // destructor
  virtual ~Matrix() { }

  // methods
  unsigned size() const { return n_; }
  virtual void multMatVec(const double * const, double * const) const = 0;       
  virtual void applyPrecond(double * const x) const { }

};

class LaplaceMatrix : public Matrix {
private:
  const unsigned size_;
public:
  // constructor
  LaplaceMatrix(const unsigned size) : Matrix(size*size), size_(size) { }

  // destructor
  ~LaplaceMatrix() { }
  unsigned operator_size() const { return size_; }
  
  // methods
  void multMatVec(const double * const, double * const) const;
};

class LaplaceMatrixJacobi : public LaplaceMatrix {
public:
  // constructor
  LaplaceMatrixJacobi(const unsigned size) 
    : LaplaceMatrix(size) { }

  // destructor
  ~LaplaceMatrixJacobi() { }

  // methods
  void applyPrecond(double * const x) const;
};  

#endif // MATRIX_HPP



