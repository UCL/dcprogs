#include <iostream>
#include <type_traits>
#include <gtest/gtest.h>

#include "../root_finder.h"
#include "random_matrix.h"
using namespace DCProgs;

// Sets up test with parameters from CH82, 1e-7 nM.
class RootFinderLowerBoundTest : public ::testing::Test {
  
  public:
  virtual void SetUp() {
    Q.resize(5, 5);
    Q <<  -3050,        50,  3000,      0,    0,
          2./3., -1502./3.,     0,    500,    0, 
             15,         0, -2065,     50, 2000, 
              0,     15000,  4000, -19000,    0, 
              0,         0,    10,      0,  -10;
  }
  protected:
    t_rmatrix Q;
};

std::ostream &operator<< (std::ostream &_stream, DCProgs::RootInterval &_interval) {
  return _stream << "@(" << _interval.start << ", "
                 << _interval.end << ") = " << _interval.multiplicity;
}

// Checks can find zero for cannonical matrix, open states.
TEST_F(RootFinderLowerBoundTest, open) {
  StateMatrix const states(Q, 2);
  DeterminantEq det(states, 1e-4, true); 
  
  EXPECT_TRUE(find_lower_bound_for_roots(det) < -3045.285776);
  EXPECT_TRUE(std::abs(find_lower_bound_for_roots(det, -4e4) + 4e4) < 1e-8);
}
TEST_F(RootFinderLowerBoundTest, close) {
  StateMatrix const states(Q, 2);
  DeterminantEq det(states, 1e-4, false); 
  
  EXPECT_TRUE(find_lower_bound_for_roots(det) < -17090.1927692368);
  EXPECT_TRUE(std::abs(find_lower_bound_for_roots(det, -2e5) + 2e5) < 1e-4);
}

class TestFindLowerBound : public ::testing::TestWithParam<t_int> {
  public:
  TestFindLowerBound() {};
};

TEST_P(TestFindLowerBound, non_singular) {

  typedef std::uniform_int_distribution<t_int> t_idist;
  
  StateMatrix Qmatrix;
  try {
    Qmatrix.matrix = nonsingular_qmatrix();
    Qmatrix.nopen = t_idist(2, Qmatrix.matrix.rows()-2)(global_mersenne());
    DeterminantEq det(Qmatrix, 1e-4, true);
    t_real const lb(  find_lower_bound_for_roots(det) );
    Eigen::EigenSolver<t_rmatrix> eigsolver(det.H(lb));
    EXPECT_TRUE((eigsolver.eigenvalues().array().real() > lb).all()) 
        << "Found eigenvalue below lower bound.\n";
  } catch(...) {
    std::cerr.precision(15);
    std::cerr << "Error for " << Qmatrix << std::endl;
    throw;
  }
}

INSTANTIATE_TEST_CASE_P(random, TestFindLowerBound, ::testing::Range(0, 300));

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

