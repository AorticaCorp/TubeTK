/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeComputeImageStatistics_h
#define __itktubeComputeImageStatistics_h

// ITK includes
#include <itkImage.h>
#include <itkImageToImageFilter.h>

namespace itk
{

namespace tube
{

/** \class ComputeImageStatistics
 * \brief Computes image statistics
 */

template< class TPixel, unsigned int VDimension >
class ComputeImageStatistics
  : public itk::ImageToImageFilter< itk::Image< float,  VDimension >,
      itk::Image< float,  VDimension > >
{
public:

  /** Standard class typedefs. */
  typedef ComputeImageStatistics                   Self;
  typedef itk::ImageToImageFilter< itk::Image< float,  VDimension >,
      itk::Image< float,  VDimension > >           Superclass;
  typedef SmartPointer< Self >                     Pointer;
  typedef SmartPointer< const Self >               ConstPointer;

  /** Custom typedefs */
  typedef itk::Image< TPixel,  VDimension >        MaskType;
  typedef itk::Image< unsigned int,  VDimension >  ConnCompType;
  typedef itk::Image< float,  VDimension >         VolumeType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ComputeImageStatistics, ImageToImageFilter );

  /** Set/Get input mask */
  itkSetObjectMacro( InputMask, MaskType );
  itkGetObjectMacro( InputMask, MaskType );

    /** Set/Get input mask */
  virtual void SetQuantiles( std::vector<float> _arg );
  itkGetMacro( Quantiles, std::vector<float> );

  /** Get Components */
  itkGetMacro( CompMean, itk::Array< double > );
  itkGetMacro( CompMin, itk::Array< double > );
  itkGetMacro( CompMax, itk::Array< double > );
  itkGetMacro( CompStdDev, itk::Array< double > );
  itkGetMacro( CompCount, itk::Array< double > );
  itkGetMacro( CompValue, itk::Array< TPixel > );
  itkGetMacro( NumberOfComponents, unsigned int );

  /** Write statistics to a CSV formatted file */
  void WriteCSVStatistics( std::string csvStatisticsFile ) const;


protected:

  ComputeImageStatistics( void );
  ~ComputeImageStatistics( void ) {};

  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Computes image statistics */
  void GenerateData( void );

private:

  typename MaskType::Pointer   m_InputMask;

  std::vector<float>           m_Quantiles;

  itk::Array< double >         m_CompMean;
  itk::Array< double >         m_CompMin;
  itk::Array< double >         m_CompMax;
  itk::Array< double >         m_CompStdDev;
  itk::Array< double >         m_CompCount;
  itk::Array< TPixel >         m_CompValue;

  unsigned int                 m_NumberOfComponents;

  vnl_matrix< double >         m_QuantileValue;

}; // End class ComputeImageStatistics

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeComputeImageStatistics.hxx"
#endif

#endif // End !defined(__itktubeComputeImageStatistics_h)
