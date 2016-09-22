/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itktubeSegmentTubesUsingMinimalPathFilter_hxx
#define __itktubeSegmentTubesUsingMinimalPathFilter_hxx

#include "itktubeSegmentTubesUsingMinimalPathFilter.h"

// MinimalPathExtraction Imports
#include "itkSpeedFunctionToPathFilter.h"
#include "itkIterateNeighborhoodOptimizer.h"
#include "itkSingleImageCostFunction.h"

#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkSelectEnclosedPoints.h"

namespace itk
{
	namespace tube
	{
		template< unsigned int Dimension, class TInputPixel > SegmentTubesUsingMinimalPathFilter< Dimension, TInputPixel >::SegmentTubesUsingMinimalPathFilter(void)
		{
			m_SpeedImage = NULL;
			m_RadiusImage = NULL;
			m_TargetTubeGroup = NULL;
			m_ConnectToTargetTubeSurface = false;
			m_NoBend = false;
			m_BendUpwards = false;
			m_BendDownwards = false;
			m_OptimizationMethod = "Regular_Step_Gradient_Descent";
			m_OptimizerTerminationValue = 2.0;
			m_OptimizerNumberOfIterations = 2000;
			m_OptimizerStepLengthFactor = 0.1;
			m_OptimizerStepLengthRelax = 0.999;
			m_StartRadius = 1;
			m_MaxRadius = 6;
			m_StepSizeForRadiusEstimation = 0.5;
			m_CostAssociatedWithExtractedTube = 0.0;
			m_Output = NULL;
		}

		template< unsigned int Dimension, class TInputPixel > void SegmentTubesUsingMinimalPathFilter< Dimension, TInputPixel >::SetIntermediatePoints(std::vector< PointType > pathPoints)
		{
			m_IntermediatePoints = pathPoints;
		} // end setintermediatepoints

		template< unsigned int Dimension, class TInputPixel > void SegmentTubesUsingMinimalPathFilter< Dimension, TInputPixel >::Update(void)
		{
			typedef itk::PolyLineParametricPath< Dimension > PathType;
			typedef itk::SpeedFunctionToPathFilter< InputImageType, PathType > PathFilterType;
			typedef itk::LinearInterpolateImageFunction< InputImageType, double >InterpolatorType;
			typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

			typedef itk::SingleImageCostFunction< InputImageType > CostFunctionType;
			typename CostFunctionType::Pointer costFunction = CostFunctionType::New();
			costFunction->SetInterpolator(interpolator);

			//Get Input image information
			typedef typename TubeType::TransformType TransformType;
			typename TransformType::InputVectorType scaleVector;
			typename TransformType::OffsetType offsetVector;
			typename InputImageType::SpacingType spacing = m_SpeedImage->GetSpacing();
			typename InputImageType::PointType origin = m_SpeedImage->GetOrigin();
			double tubeSpacing[Dimension];

			for (unsigned int k = 0; k < Dimension; ++k)
			{
				scaleVector[k] = spacing[k];
				offsetVector[k] = origin[k];
				tubeSpacing[k] = spacing[k];
			}
			// Create path filter & path information
			typename PathFilterType::Pointer pathFilter = PathFilterType::New();
			pathFilter->SetInput(m_SpeedImage);
			pathFilter->SetCostFunction(costFunction);
			pathFilter->SetTerminationValue(m_OptimizerTerminationValue);

			typedef itk::SpeedFunctionPathInformation< PointType > PathInformationType;
			typename PathInformationType::Pointer pathInfo = PathInformationType::New();
			pathInfo->SetStartPoint(m_StartPoint);
			for (unsigned int i = 0; i < m_IntermediatePoints.size(); i++)
			{
				pathInfo->AddWayPoint(m_IntermediatePoints[i]);
			}
			double radiusAtMin = 0.0;

			// if user provides an end point for centerline
			if (m_EndPoint.Size())
			{
				pathInfo->SetEndPoint(m_EndPoint);
			} //end check - if user provides an end point for centerline

			// if user does not provide an end point but instead provides a target tube and says there is no bend in vessel
			// in this case use the smart closest point search that @vicory implemented
			// (smart closest point search - user does not provide an end point for centerline and does not say there is a bend in vessel.
			//  so centerline has to find an end point that is inside the target tube in a direction as straight as possible. this is the search that @vicory implemented)
			if (m_TargetTubeGroup && m_NoBend)
			{
				std::cout << "no bend in vessel but i see a target tube therefore have to clip" << std::endl;

				// Find point on target tube w/ earliest arrival and use that as end point
				typedef FastMarchingUpwindGradientImageFilter< InputImageType, InputImageType > FastMarchingType;

				typedef typename FastMarchingType::NodeContainer NodeContainer;
				typedef typename FastMarchingType::NodeType NodeType;
				typename FastMarchingType::Pointer marching = FastMarchingType::New();
				marching->SetInput(m_SpeedImage);
				marching->SetGenerateGradientImage(false);
				marching->SetTargetOffset(2.0 * m_OptimizerTerminationValue);
				marching->SetTargetReachedModeToAllTargets();

				// Use closest point in Euclidean sense as initial guess
				PointType pointPath;
				this->IsPointTooNear(m_TargetTubeGroup, m_StartPoint, pointPath);

				typename PathInformationType::Pointer tempPath = PathInformationType::New();
				tempPath->SetStartPoint(m_StartPoint);
				tempPath->SetEndPoint(pointPath);

				typename InputImageType::IndexType indexTargetPrevious, indexTargetNext;
				m_SpeedImage->TransformPhysicalPointToIndex(tempPath->PeekPreviousFront(), indexTargetPrevious);
				m_SpeedImage->TransformPhysicalPointToIndex(tempPath->PeekNextFront(), indexTargetNext);

				NodeType nodeTargetPrevious;
				NodeType nodeTargetNext;
				nodeTargetPrevious.SetValue(0.0);
				nodeTargetNext.SetValue(0.0);
				nodeTargetPrevious.SetIndex(indexTargetPrevious);
				nodeTargetNext.SetIndex(indexTargetNext);
				typename NodeContainer::Pointer targets = NodeContainer::New();
				targets->Initialize();
				targets->InsertElement(0, nodeTargetPrevious);
				targets->InsertElement(1, nodeTargetNext);
				marching->SetTargetPoints(targets);

				typename InputImageType::IndexType indexTrial;
				m_SpeedImage->TransformPhysicalPointToIndex(tempPath->GetCurrentFrontAndAdvance(), indexTrial);
				NodeType nodeTrial;
				nodeTrial.SetValue(0.0);
				nodeTrial.SetIndex(indexTrial);
				typename NodeContainer::Pointer trial = NodeContainer::New();
				trial->Initialize();
				trial->InsertElement(0, nodeTrial);
				marching->SetTrialPoints(trial);

				marching->UpdateLargestPossibleRegion();
				typename InputImageType::Pointer arrival = marching->GetOutput();

				typename InputSpatialObjectType::ChildrenListPointer tubeList = m_TargetTubeGroup->GetChildren();
				typename TubeType::Pointer curTube = dynamic_cast<TubeType*>(tubeList->begin()->GetPointer());
				typename TubeType::PointListType pointList = curTube->GetPoints();
				curTube->ComputeObjectToWorldTransform();

				double minDist = itk::NumericTraits<double>::max();
				int minInd = -1;
				for (int i = 0; i < pointList.size(); i++)
				{
					typename InputImageType::IndexType ind;
					ind[0] = pointList[i].GetPosition()[0];
					ind[1] = pointList[i].GetPosition()[1];
					ind[2] = pointList[i].GetPosition()[2];

					double d = arrival->GetPixel(ind);
					if (d < minDist)
					{
						minDist = d;
						minInd = i;
						radiusAtMin = pointList[i].GetRadius();
					}
				}

				pathInfo->SetEndPoint(curTube->GetIndexToWorldTransform()->TransformPoint(pointList[minInd].GetPosition()));
			} //end check - if user does not provide an end point but instead provides a target tube and says there is no bend in vessel

			// if user does not provide an end point but provides a target tube and says that the branch vessel bends upwards
			// in this case set the end point for the vessel's centerkine to the top  of the target tube
			if (m_TargetTubeGroup && m_BendUpwards)
			{
				std::cout << "bend upwards. therefore end point for centerline is top pt of target tube and i see a target tube therefore have to clip" << std::endl;
				typename InputSpatialObjectType::ChildrenListPointer tubeList = m_TargetTubeGroup->GetChildren();
				typename TubeType::Pointer curTube = dynamic_cast<TubeType*>(tubeList->begin()->GetPointer());
				typename TubeType::PointListType pointList = curTube->GetPoints();
				curTube->ComputeObjectToWorldTransform();

				pathInfo->SetEndPoint(curTube->GetIndexToWorldTransform()->TransformPoint(pointList[pointList.size() - 1].GetPosition()));

				PointType tmpPt = curTube->GetIndexToWorldTransform()->TransformPoint(pointList[pointList.size() - 1].GetPosition());
				std::cout << tmpPt[0] << " " << tmpPt[1] << " " << tmpPt[2] << std::endl;
			} // end check - if user does not provide an end point but provides a target tube and says that the branch vessel bends upward

			if (m_TargetTubeGroup && m_BendDownwards)
			{
				std::cout << "bend downwards. therefore end point for centerline is bottom pt of target tube and i see a target tube therefore have to clip" << std::endl;
				typename InputSpatialObjectType::ChildrenListPointer tubeList = m_TargetTubeGroup->GetChildren();
				typename TubeType::Pointer curTube = dynamic_cast<TubeType*>(tubeList->begin()->GetPointer());
				typename TubeType::PointListType pointList = curTube->GetPoints();
				curTube->ComputeObjectToWorldTransform();

				pathInfo->SetEndPoint(curTube->GetIndexToWorldTransform()->TransformPoint(pointList[0].GetPosition()));

				PointType tmpPt = curTube->GetIndexToWorldTransform()->TransformPoint(pointList[0].GetPosition());
				std::cout << tmpPt[0] << " " << tmpPt[1] << " " << tmpPt[2] << std::endl;

			} // end check - if user does not provide an end point but provides a target tube and says that the branch vessel bends upward

			pathFilter->AddPathInformation(pathInfo);

			// Set Optimizer
			if (m_OptimizationMethod == "Iterate_Neighborhood")
			{
				// Create IterateNeighborhoodOptimizer
				typedef itk::IterateNeighborhoodOptimizer OptimizerType;
				typename OptimizerType::Pointer optimizer = OptimizerType::New();
				optimizer->MinimizeOn();
				optimizer->FullyConnectedOn();
				typename OptimizerType::NeighborhoodSizeType size(Dimension);
				for (unsigned int i = 0; i < Dimension; i++)
				{
					size[i] = m_SpeedImage->GetSpacing()[i] * m_OptimizerStepLengthFactor;
				}
				optimizer->SetNeighborhoodSize(size);
				pathFilter->SetOptimizer(optimizer);
			}
			else if (m_OptimizationMethod == "Gradient_Descent")
			{
				// Create GradientDescentOptimizer
				typedef itk::GradientDescentOptimizer OptimizerType;
				typename OptimizerType::Pointer optimizer = OptimizerType::New();
				optimizer->SetNumberOfIterations(m_OptimizerNumberOfIterations);
				pathFilter->SetOptimizer(optimizer);
			}
			else if (m_OptimizationMethod == "Regular_Step_Gradient_Descent")
			{
				// Compute the minimum spacing
				double minspacing = spacing[0];
				for (unsigned int dim = 0; dim < Dimension; dim++)
				{
					if (spacing[dim] < minspacing)
					{
						minspacing = spacing[dim];
					}
				}
				// Create RegularStepGradientDescentOptimizer
				typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
				typename OptimizerType::Pointer optimizer = OptimizerType::New();
				optimizer->SetNumberOfIterations(m_OptimizerNumberOfIterations);
				optimizer->SetMaximumStepLength(1.0 * m_OptimizerStepLengthFactor * minspacing);
				optimizer->SetMinimumStepLength(0.5 * m_OptimizerStepLengthFactor * minspacing);
				optimizer->SetRelaxationFactor(m_OptimizerStepLengthRelax);
				pathFilter->SetOptimizer(optimizer);
			}
			else
			{
				return;
			}

			try
			{
				pathFilter->Update();
			}
			catch (itk::ExceptionObject & err)
			{
				std::stringstream out;
				out << "ExceptionObject caught !" << std::endl;
				out << err << std::endl;
				return;
			}

			// Create output TRE file
			m_Output = InputSpatialObjectType::New();

			// Update tubes transform
			m_Output->GetObjectToParentTransform()->SetScale(scaleVector);
			m_Output->GetObjectToParentTransform()->SetOffset(offsetVector);
			m_Output->GetObjectToParentTransform()->SetMatrix(m_SpeedImage->GetDirection());
			m_Output->ComputeObjectToWorldTransform();
			m_CostAssociatedWithExtractedTube = 0.0;
			for (unsigned int i = 0; i < pathFilter->GetNumberOfOutputs(); i++)
			{
				// Get the path
				typename PathType::Pointer path = pathFilter->GetOutput(i);
				// Check path is valid
				if (path->GetVertexList()->Size() == 0)
				{
					std::cout << "WARNING: Path " << (i + 1) << " contains no points!" << std::endl;
					continue;
				}

				// Output centerline in TRE file
				typename TubeType::PointListType tubePointList;
				typename PathType::VertexListType * vertexList = path->GetVertexList();
				unsigned int k, ctr = 0;

				for (k = 0; k < vertexList->Size(); k++)
				{
					PointType pathPoint;
					m_SpeedImage->TransformContinuousIndexToPhysicalPoint(vertexList->GetElement(k), pathPoint);

					typename InputImageType::IndexType imageIndex;
					m_SpeedImage->TransformPhysicalPointToIndex(pathPoint, imageIndex);

					if (m_SpeedImage->TransformPhysicalPointToIndex(pathPoint, imageIndex))
					{
						m_CostAssociatedWithExtractedTube += m_SpeedImage->GetPixel(imageIndex);
					}

					TubePointType tubePoint;
					tubePoint.SetPosition(vertexList->GetElement(k));
					tubePoint.SetRadius(m_SpeedImage->GetPixel(imageIndex));
					tubePoint.SetID(k);
					tubePointList.push_back(tubePoint);
				} // end loop: number of points in current tube

				if (m_ConnectToTargetTubeSurface)
				{
					typename TubeType::PointListType tubePointList;

					for (int ii = 0; ii < vertexList->Size(); ii++) // TODO: replace this loop with loop through tubepointlist 
					{
						PointType pathPoint;
						m_SpeedImage->TransformContinuousIndexToPhysicalPoint(vertexList->GetElement(ii), pathPoint);

						typename InputImageType::IndexType imageIndex;
						m_SpeedImage->TransformPhysicalPointToIndex(pathPoint, imageIndex);

						PointType nearPoint;
						//bool isnear = this->IsPointTooNear(m_TargetTubeGroup, pathPoint, nearPoint);

						typename InputSpatialObjectType::ChildrenListPointer tubeList = m_TargetTubeGroup->GetChildren();
						typename TubeType::Pointer curTube = dynamic_cast<TubeType*>(tubeList->begin()->GetPointer());
						typename TubeType::PointListType pointList = curTube->GetPoints();
						curTube->ComputeObjectToWorldTransform();

						//if (!curTube->IsInside(pathPoint))
						if(!this->IsPointTooNear(m_TargetTubeGroup, pathPoint, nearPoint))
						{	
							TubePointType tubePoint;
							tubePoint.SetPosition(vertexList->GetElement(ii));
							tubePoint.SetRadius(m_SpeedImage->GetPixel(imageIndex));
							tubePoint.SetID(ctr);
							tubePointList.push_back(tubePoint);
							ctr++;
						}
					}
					typename TubeType::Pointer pTube = TubeType::New();
					pTube->SetPoints(tubePointList);
					pTube->ComputeTangentAndNormals();
					pTube->SetSpacing(tubeSpacing);
					pTube->SetId(i);

					m_Output->AddSpatialObject(pTube);
					m_Output->ComputeObjectToWorldTransform();
				} // end check: user wants to clip branch with target tube (note for SPMT, target tube = aorta)
				else
				{
					typename TubeType::Pointer pTube = TubeType::New();
					pTube->SetPoints(tubePointList);
					pTube->ComputeTangentAndNormals();
					pTube->SetSpacing(tubeSpacing);
					pTube->SetId(i);

					m_Output->AddSpatialObject(pTube);
					m_Output->ComputeObjectToWorldTransform();
				}

			} // end loop: number of outputs from fast marching filter (note: number of outputs=1 for SPMT)


			//// Extract Radius
			//if( m_RadiusImage )
			//  {
			//  typedef itk::tube::RadiusExtractor2< InputImageType > RadiusExtractorType;
			//  typename RadiusExtractorType::Pointer radiusExtractor = RadiusExtractorType::New();
			//  radiusExtractor->SetInputImage( m_RadiusImage );
			//  radiusExtractor->SetRadiusStart( m_StartRadius );
			//  radiusExtractor->SetRadiusMin( 0.2 );
			//  radiusExtractor->SetRadiusMax( m_MaxRadius );
			//  radiusExtractor->SetRadiusStep( m_StepSizeForRadiusEstimation );
			//  radiusExtractor->SetRadiusTolerance( 0.025 );
			//  radiusExtractor->SetDebug( false );
			//  radiusExtractor->ExtractRadii( pTube );
			//  }


		} // end update

		template< unsigned int Dimension, class TInputPixel > bool SegmentTubesUsingMinimalPathFilter< Dimension, TInputPixel >::IsPointTooNear(const InputSpatialObjectType * sourceTubeGroup, PointType outsidePoint, PointType &nearestPoint)
		{

			double minDistance = itk::NumericTraits<double>::max();
			double nearestPointRadius = 0.0;

			typename InputSpatialObjectType::ChildrenListPointer sourceTubeList = sourceTubeGroup->GetChildren();
			for (typename InputSpatialObjectType::ChildrenListType::iterator tubeList_it = sourceTubeList->begin(); tubeList_it != sourceTubeList->end(); ++tubeList_it)
			{
				//**** Source Tube **** :
				typename TubeType::Pointer pCurSourceTube = dynamic_cast<TubeType*>(tubeList_it->GetPointer());
				//dynamic_cast verification
				if (!pCurSourceTube)
				{
					return EXIT_FAILURE;
				}
				pCurSourceTube->ComputeObjectToWorldTransform();
				//Get points in current source tube
				typename TubeType::PointListType pointList = pCurSourceTube->GetPoints();
				//Get Index to World Transformation
				typename TubeType::TransformType * pTubeIndexPhysTransform = pCurSourceTube->GetIndexToWorldTransform();

				for (typename TubeType::PointListType::const_iterator pointList_it = pointList.begin(); pointList_it != pointList.end(); ++pointList_it)
				{
					TubePointType curSourcePoint = *pointList_it;
					//Transform parameters in physical space
					typename TubePointType::PointType curSourcePos = pTubeIndexPhysTransform->TransformPoint(curSourcePoint.GetPosition());
					double distance = curSourcePos.SquaredEuclideanDistanceTo(outsidePoint);

					if (minDistance > distance)
					{
						minDistance = distance;
						for (unsigned int i = 0; i < Dimension; i++)
						{
							nearestPoint[i] = curSourcePos[i];
						}
						nearestPointRadius = curSourcePoint.GetRadius();
					}
				}
			}

			//if ( minDistance < nearestPointRadius*nearestPointRadius)
			if (sqrt(minDistance) - nearestPointRadius < -1e-4) 
			{
				return true;
			}
			else
			{
				return false;
			}
		}

		template< unsigned int Dimension, class TInputPixel > void SegmentTubesUsingMinimalPathFilter< Dimension, TInputPixel >::PrintSelf(std::ostream & os, Indent indent) const
		{
			Superclass::PrintSelf(os, indent);

		}

	} // end namespace tube
} // end namespace itk

#endif // End !defined( __itktubeSegmentTubesUsingMinimalPathFilter_hxx )
