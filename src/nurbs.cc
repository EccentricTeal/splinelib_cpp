#include "splinelib_cpp/nurbs.hh"

namespace splinelib
{
  NurbsData::NurbsData( uint16_t dimension ):
  nurbs_dimension_( dimension )
  {
    ;
  }


  NurbsData::~NurbsData()
  {
    ;
  }


  void NurbsData::setControlPoints( std::vector<NurbsControlPoint>& new_cp )
  {
    control_points_.clear();
    for( auto itr = new_cp.begin(); itr != new_cp.end(); ++itr )
    { control_points_.push_back( *itr ); }
  }


  void NurbsData::setKnotVector( std::vector<double>& new_kv )
  {
    knot_vector_.clear();
    for( auto itr = new_kv.begin(); itr != new_kv.end(); ++itr )
    { knot_vector_.push_back( *itr ); }
  }


  NurbsError NurbsData::evalCurveData( void )
  {
    //empty evaluation
    if( control_points_.empty() == true )
    { return NurbsError::Err_ControlPointEmpty; }
    if( knot_vector_.empty() == true )
    { return NurbsError::Err_KnotVectorEmpty; }

    //data number evaluation
    if( knot_vector_.size() < 3 )
    { return NurbsError::Err_KnotVectorInvalid; }

    for( auto itr = std::next( knot_vector_.begin() ); itr != knot_vector_.end(); ++itr )
    {
      if( *itr < *( std::prev( itr ) ) ){ return NurbsError::Err_KnotVectorInvalid; }
    }
    
    return NurbsError::NoError;
  }


  Eigen::Vector3d NurbsData::getPointOnCurve( double t )
  {
    double l = getKnotSectionBeginNum( t );
    uint16_t n = nurbs_dimension_;

    double w = calcW( t, l, n );
    Eigen::Vector3d wP = calcwP( t, l, n );

    Eigen::Vector3d P;
    if( w != 0.0 ){ P = wP / w; }

    return wP;
  }


  Eigen::Vector3d NurbsData::getDifferential( double t )
  {
    double l = getKnotSectionBeginNum( t );
    uint16_t n = nurbs_dimension_;

    double w1 = calcW( t, l-1, n-1 );
    double w2 = calcW( t, l, n-1 );
    double w3 = calcW( t, l, n );
    Eigen::Vector3d wP1 = calcwP( t, l-1, n-1 );
    Eigen::Vector3d wP2 = calcwP( t, l, n-1 );

    Eigen::Vector3d diff;
    if( w3 != 0.0 )
    {
      diff = 
        ( n * w1 * w2 / ( w3 * w3 ) ) *
        ( ( ( wP1 / w1 ) - ( wP2 / w2 ) ) / ( knot_vector_[l+1] - knot_vector_[l] ) );
    }

    return diff;
  }


  uint16_t NurbsData::getKnotSectionBeginNum( double t )
  {
    uint16_t vector_position = 0;

    if( *( knot_vector_.end()-1 ) == t ){ vector_position =  knot_vector_.size(); }
    else if( *( knot_vector_.begin() ) == t )
    {
      while( knot_vector_[ vector_position+1 ] == t && (vector_position+1) < knot_vector_.size() )
      {
        vector_position++;
      }
    }
    else
    {
      for( auto itr = knot_vector_.begin()+1; itr != knot_vector_.end(); ++itr )
      {
        if( *itr > t ){ break; }
        vector_position++;
      }
    }

    return vector_position;
  }


  double NurbsData::calcAlpha( double t, uint16_t i, uint16_t k )
  {
    uint16_t n = nurbs_dimension_;
    double alpha = 0.0;

    if( ( knot_vector_[( i + n + 1 - k )] - knot_vector_[i] ) != 0.0 )
    {
      alpha = ( t - knot_vector_[i] ) / ( knot_vector_[( i + n + 1 - k )] - knot_vector_[i] );
    }
    
    return alpha;
  }


  double NurbsData::calcW( double t, double i, double k )
  {
    double w = 0.0;

    if( k == 0 )
    { w = control_points_[i].weight; }
    else
    {
      w = ( 1.0 - calcAlpha( t, i, k ) ) * calcW( t, i-1, k-1 ) +
          calcAlpha( t, i, k ) * calcW( t, i, k-1 );
    }

    return w;
  }


  Eigen::Vector3d NurbsData::calcwP( double t, double i, double k )
  {
    Eigen::Vector3d wP;

    if( k == 0)
    {
      wP = control_points_[i].weight * control_points_[i].point;
    }
    else
    {
      double alpha_1 = 1.0 - calcAlpha( t, i, k );
      Eigen::Vector3d wP1 = calcwP( t, i-1, k-1 );

      double alpha_2 = calcAlpha( t, i, k );
      Eigen::Vector3d wP2 = calcwP( t, i, k-1 );

      wP = ( alpha_1 * wP1 ) + ( alpha_2 * wP2 );
    }
    
    return wP;
  }

}