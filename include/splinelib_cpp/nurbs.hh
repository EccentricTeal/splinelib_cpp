#ifndef SPLINELIB_CPP__NURBS__HH
#define SPLINELIB_CPP__NURBS__HH

#include "Eigen/Dense"
#include <vector>

namespace splinelib
{
  struct NurbsControlPoint
  {
    Eigen::Vector3d point;
    double weight = 1.0;
  };

  enum class NurbsError
  {
    NoError = 0,
    Err_ControlPointEmpty = -1,
    Err_KnotVectorEmpty = -2,
    Err_KnotVectorInvalid = -3
  };

  /* Mathematical theory is included from a website next:
 * https://www.f-sp.com/entry/2016/08/07/020613
 */
  class NurbsData
  {
    /** Constants **/
    private:

    /** Private Objects **/
    //NURBS Data
    private:
      const uint16_t nurbs_dimension_;
      std::vector<NurbsControlPoint> control_points_;
      std::vector<double> knot_vector_;

    /** Constructor, Destructor **/
    public:
      NurbsData( uint16_t dimension );
      ~NurbsData();

    /** Public Methods **/
    public:
      void setControlPoints( std::vector<NurbsControlPoint>& new_cp );
      void setKnotVector( std::vector<double>& new_kv );
      NurbsError evalCurveData( void );
      Eigen::Vector3d getPointOnCurve( double t );
      Eigen::Vector3d getDifferential( double t );
    
    /** Private Methods **/
    private:
      uint16_t getKnotSectionBeginNum( double t );
      double calcAlpha( double t, uint16_t i, uint16_t k );
      double calcW( double t, double i, double k  );
      Eigen::Vector3d calcwP( double t, double i, double k );
  };
}



#endif