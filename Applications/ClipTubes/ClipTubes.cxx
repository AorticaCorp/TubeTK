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
#include <iostream>
#include <sstream>

#include "itkTimeProbesCollectorBase.h"
#include "tubeMessage.h"
#include "tubeMacro.h"
#include "metaScene.h"

#include "itkGroupSpatialObject.h"
#include "itkImageFileReader.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkTimeProbesCollectorBase.h"
#include "tubeCropTubes.h"
#include "ClipTubesCLP.h"

template< unsigned int DimensionT >
int DoIt (int argc, char * argv[])
{
  PARSE_ARGS;

  // Ensure that the input image dimension is valid
  // We only support 2D and 3D Images due to the
  // limitation of itkTubeSpatialObject
  if( DimensionT != 2 && DimensionT != 3 )
    {
    tube::ErrorMessage(
      "Error: Only 2D and 3D data is currently supported.");
    return EXIT_FAILURE;
    }

  // The timeCollector to perform basic profiling of algorithmic components
  itk::TimeProbesCollectorBase timeCollector;

  // Load TRE File
  tubeStandardOutputMacro( << "\n>> Loading TRE File" );

  typedef itk::SpatialObjectReader< DimensionT >          TubesReaderType;
  typedef double                                          PixelType;
  typedef itk::Vector< double, DimensionT >               VectorType;
  typedef itk::Point< double, DimensionT >                PointType;

  timeCollector.Start( "Loading Input TRE File" );

  typedef tube::CropTubes< DimensionT > FilterType;
  typename FilterType::Pointer filter = FilterType::New();


  typename TubesReaderType::Pointer tubeFileReader = TubesReaderType::New();

  try
    {
    tubeFileReader->SetFileName( inputTREFile.c_str() );
    tubeFileReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error loading TRE File: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  filter->SetInput( tubeFileReader->GetGroup() );

  timeCollector.Stop( "Loading Input TRE File" );

  timeCollector.Start( "Loading Input Parameters" );


  timeCollector.Stop( "Writing output TRE file" );

  timeCollector.Report();
  return EXIT_SUCCESS;
}


// Main
int main( int argc, char * argv[] )
{
  try
    {
    PARSE_ARGS;
    }
  catch( const std::exception & err )
    {
    tube::ErrorMessage( err.what() );
    return EXIT_FAILURE;
    }
  PARSE_ARGS;
  
  MetaScene *mScene = new MetaScene;
  mScene->Read( inputTREFile.c_str() );

  if( mScene->GetObjectList()->empty() )
    {
    tubeWarningMacro( << "Input TRE file has no spatial objects" );
    delete mScene;
    return EXIT_SUCCESS;
    }

  switch( mScene->GetObjectList()->front()->NDims() )
    {
    case 2:
      {
      bool result = DoIt<2>( argc, argv );
      delete mScene;
      return result;
      break;
      }
    case 3:
      {
      bool result = DoIt<3>( argc, argv );
      delete mScene;
      return result;
      break;
      }
    default:
      {
      tubeErrorMacro(
        << "Error: Only 2D and 3D data is currently supported." );
      delete mScene;
      return EXIT_FAILURE;
      break;
      }
    }
  return EXIT_FAILURE;
}
