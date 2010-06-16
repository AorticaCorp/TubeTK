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
#ifndef __itkJointHistogramImageFunction_txx
#define __itkJointHistogramImageFunction_txx

#include "itkMinimumMaximumImageCalculator.h"
#include "itkAddImageFilter.h"
#include "itkSquareImageFilter.h"

#include "itkJointHistogramImageFunction.h"


namespace itk
{

/**
 * Set the input Image
 */
template <class TInputImage, class TCoordRep>
JointHistogramImageFunction<TInputImage,TCoordRep>
::JointHistogramImageFunction()
{
  m_InputMask = 0;
  m_FeatureWidth = 20;
  m_SumHistogram = 0;
  m_SumOfSquaresHistogram = 0;
  m_NumberOfSamples = 0;
  m_NumberOfComputedSamples = 0;
  m_ImageMin = 0;
  m_ImageMax = 0;
  m_MaskMin = 0;
  m_MaskMax = 0;
  this->SetHistogramSize( 20 );
}

template <class TInputImage, class TCoordRep>
void 
JointHistogramImageFunction<TInputImage,TCoordRep>
::SetHistogramSize( const unsigned int& size )
{
  m_HistogramSize = size;
  typename HistogramType::IndexType start;
  typename HistogramType::SizeType histSize;
  start[0] = start[1] = 0;
  histSize[0] = histSize[1] = size;
  typename HistogramType::RegionType region;
  region.SetIndex( start );
  region.SetSize( histSize );

  m_SumHistogram = HistogramType::New();
  m_SumHistogram->SetRegions( region );
  m_SumHistogram->Allocate();
  m_SumHistogram->FillBuffer( 0 );

  m_SumOfSquaresHistogram = HistogramType::New();
  m_SumOfSquaresHistogram->SetRegions( region );
  m_SumOfSquaresHistogram->Allocate();
  m_SumOfSquaresHistogram->FillBuffer( 0 );

  m_NumberOfComputedSamples = m_NumberOfSamples = 0;
}

template <class TInputImage, class TCoordRep>
void
JointHistogramImageFunction<TInputImage,TCoordRep>
::SetInputImage( const InputImageType * ptr )
{
  this->Superclass::SetInputImage( ptr );

  typedef itk::MinimumMaximumImageCalculator<InputImageType> MinMaxType;
  typename MinMaxType::Pointer calculator = MinMaxType::New();
  calculator->SetImage( this->GetInputImage() );
  calculator->Compute();
  m_ImageMin = calculator->GetMinimum();
  m_ImageMax = calculator->GetMaximum();
}

template <class TInputImage, class TCoordRep>
void
JointHistogramImageFunction<TInputImage,TCoordRep>
::SetInputMask( const typename InputImageType::Pointer mask )
{
  m_InputMask = mask;

  typedef itk::MinimumMaximumImageCalculator<InputImageType> MinMaxType;
  typename MinMaxType::Pointer calculator = MinMaxType::New();
  calculator->SetImage( m_InputMask );
  calculator->Compute();
  m_MaskMin = calculator->GetMinimum();
  m_MaskMax = calculator->GetMaximum();
}

template <class TInputImage, class TCoordRep>
double 
JointHistogramImageFunction<TInputImage,TCoordRep>
::EvaluateAtIndex( const IndexType & index ) const
{
  if( m_NumberOfComputedSamples < m_NumberOfSamples )
    {
    
    m_NumberOfComputedSamples = m_NumberOfSamples;
    }
  PixelType value = this->GetInputImage()->GetPixel( index );
  return static_cast<double>( value );
}

template <class TInputImage, class TCoordRep>
void 
JointHistogramImageFunction<TInputImage,TCoordRep>
::PrecomputeAtIndex( const IndexType & index )
{
  typename HistogramType::Pointer hist;
  this->ComputeHistogramAtIndex( index, hist );

  typedef itk::AddImageFilter< HistogramType, HistogramType, HistogramType>
                                                             AdderType;
  typedef itk::SquareImageFilter< HistogramType, HistogramType >                                          
                                                             SquareType;  

  typename AdderType::Pointer adder = AdderType::New();
  typename SquareType::Pointer square = SquareType::New();
      
  // Calculate Running Sum
  adder->SetInput1( m_SumHistogram );
  adder->SetInput2( hist );
  adder->Update();
  m_SumHistogram = adder->GetOutput();
  
  // Calculate Running Sum of Squares
  adder = AdderType::New();
  square->SetInput( hist );
  square->Update();
  adder->SetInput1( m_SumOfSquaresHistogram );
  adder->SetInput2( square->GetOutput() );
  adder->Update();
  m_SumOfSquaresHistogram = adder->GetOutput();
  
  ++m_NumberOfSamples;
}

template <class TInputImage, class TCoordRep>
void 
JointHistogramImageFunction<TInputImage,TCoordRep>
::ComputeMeanAndStandardDeviation()
{
  
}

template <class TInputImage, class TCoordRep>
void
JointHistogramImageFunction<TInputImage,TCoordRep>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  os << indent << "m_InputMask = " << m_InputMask << std::endl;
  os << indent << "m_SumHistogram = " << m_SumHistogram << std::endl;
  os << indent << "m_SumOfSquaresHistogram = "
     << m_SumOfSquaresHistogram << std::endl;
  os << indent << "m_MeanHistogram = " << m_MeanHistogram << std::endl;
  os << indent << "m_StandardDeviationHistogram = "
     << m_StandardDeviationHistogram << std::endl;
  os << indent << "m_FeatureWidth = " << m_FeatureWidth << std::endl;
  os << indent << "m_HistogramSize = " << m_HistogramSize << std::endl;
  os << indent << "m_NumberOfSamples = " << m_NumberOfSamples << std::endl;
  os << indent << "m_NumberOfComputedSamples = "
     << m_NumberOfComputedSamples << std::endl;
  os << indent << "m_ImageMin = " << m_InputMask << std::endl;
  os << indent << "m_ImageMax = " << m_InputMask << std::endl;
  os << indent << "m_MaskMin = " << m_MaskMin << std::endl;
  os << indent << "m_MaskMax = " << m_MaskMax << std::endl;
}

template <class TInputImage, class TCoordRep>
void
JointHistogramImageFunction<TInputImage,TCoordRep>
::ComputeHistogramAtIndex( const IndexType& index, typename HistogramType::Pointer hist ) const
{
  typename InputImageType::PixelType minInput = m_ImageMin;
  typename InputImageType::PixelType maxInput = m_ImageMax;
  typename InputImageType::PixelType normInput = 0 - minInput;
  typename InputImageType::PixelType rangeInput = maxInput - minInput;
  typename InputImageType::PixelType stepInput = rangeInput / m_HistogramSize;
  typename InputImageType::PixelType minMask = m_MaskMin;
  typename InputImageType::PixelType maxMask = m_MaskMax;
  typename InputImageType::PixelType normMask = 0 - minMask;
  typename InputImageType::PixelType rangeMask = maxMask - minMask;
  typename InputImageType::PixelType stepMask = rangeMask / m_HistogramSize;

  typename HistogramType::IndexType histStart;
  typename HistogramType::SizeType histSize;
  histStart[0] = histStart[1] = 0;
  histSize[0] = histSize[1] = m_HistogramSize;
  typename HistogramType::RegionType histRegion;
  histRegion.SetIndex( histStart );
  histRegion.SetSize( histSize );
  hist = HistogramType::New();
  hist->SetRegions( histRegion );
  hist->Allocate();
  hist->FillBuffer( 0 );

  typedef itk::ImageRegionConstIterator<InputImageType> ConstIteratorType;
  typedef itk::ImageRegionIterator<InputImageType>      IteratorType;

  typename InputImageType::SizeType size;
  typename InputImageType::RegionType region;
  IndexType origin;
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    origin[i] = index[i] - ( m_FeatureWidth * 0.5 );
    size[i] = m_FeatureWidth;
    }
  region.SetSize( size );
  region.SetIndex( origin );

  ConstIteratorType inputItr( this->GetInputImage(), region );
  ConstIteratorType maskItr( m_InputMask, region );
  while( !inputItr.IsAtEnd() && !maskItr.IsAtEnd() )
    {
    typename HistogramType::IndexType cur;
    cur[0] = ( ( inputItr.Get() + normInput ) / stepInput );
    cur[1] = ( ( maskItr.Get() + normMask ) / stepMask );
    if( cur[0] > (int)(m_HistogramSize) - 1 )
      {
      cur[0] = (int)(m_HistogramSize) - 1;
      }
    if( cur[1] > (int)(m_HistogramSize) - 1 )
      {
      cur[1] = (int)(m_HistogramSize) - 1;
      }
    if( cur[0] < 0 )
      {
      cur[0] = 0;
      }
    if( cur[1] < 0 )
      {
      cur[1] = 0;
      }

    hist->SetPixel(cur, hist->GetPixel(cur) + 1);
    
    ++inputItr;
    ++maskItr;
    }

  this->GetInputImage()->GetPixel( index );
}

} // end namespace itk

#endif