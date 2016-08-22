//
//  PolyDataReceiver.cpp
//  Hololens polydata receiver
//
//  Created by Longquan Chen on 8/22/16.
//
//

#include <stdio.h>

#include <iostream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <cstring>
//OpenIGTLink Include
#include "igtlOSUtil.h"
#include "igtlClientSocket.h"
#include "igtlMessageHeader.h"
#include "igtlPolyDataMessage.h"
#include "igtlTransformMessage.h"
#include "igtlMultiThreader.h"
#include "igtlConditionVariable.h"
#include "igtlTimeStamp.h"

// VTK includes
#include <vtkNew.h>
#include <vtkCallbackCommand.h>
#include <vtkImageData.h>
#include <vtkTransform.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPolyLine.h>
#include <vtkPolygon.h>
#include <vtkVertex.h>
#include <vtkTriangleStrip.h>
#include <vtkFloatArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkFieldData.h>


#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkVector.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkInteractorStyle.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>

vtkSmartPointer<vtkUnsignedCharArray> colors;
vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter;
vtkSmartPointer<vtkPoints>  cloud;
vtkSmartPointer<vtkPolyData> polyData;
vtkSmartPointer<vtkPolyDataMapper> mapper;
vtkSmartPointer<vtkActor> actor;
vtkSmartPointer<vtkRenderer> renderer;
vtkSmartPointer<vtkRenderWindow> renderWindow;
vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
bool interactionActive;
igtl::ConditionVariable::Pointer conditionVar;
igtl::SimpleMutexLock * localMutex;
int frameNum = 0;

bool ReceivePolyDataStream(igtl::Socket * socket, igtl::MessageHeader::Pointer header)
{
  //std::cerr << "Receiving TRANSFORM data type. " << header->GetDeviceName()<< std::endl;
  
  // Create a message buffer to receive transform data
  igtl::MessageBase::Pointer buffer = igtl::MessageBase::New();
  buffer->SetMessageHeader(header);
  buffer->AllocatePack();
  // Receive transform data from the socket
  //igtl::TimeStamp::Pointer timer = igtl::TimeStamp::New();
  //long timePre = timer->GetTimeStampInNanoseconds();
  int read = socket->Receive(buffer->GetPackBodyPointer(), buffer->GetPackBodySize());
  std::cerr<< " read " << buffer->GetPackBodySize() <<std::endl;
  if (read != buffer->GetPackBodySize())
  {
    std::cerr<< "Only read " << read << " but expected to read "
                   << buffer->GetPackBodySize() << "\n";
  }
  // Create a message buffer to receive image data
  igtl::PolyDataMessage::Pointer polyDataMsg;
  polyDataMsg = igtl::PolyDataMessage::New();
  polyDataMsg->Copy(buffer); // !! TODO: copy makes performance issue.
  
  // Deserialize the data
  // If CheckCRC==0, CRC check is skipped.
  int c = polyDataMsg->Unpack(0);
  
  if ((c & igtl::MessageHeader::UNPACK_BODY) == 0) // if CRC check fails
  {
    std::cerr<<"Unable to create MRML node from incoming POLYDATA message. Failed to unpack the message";
    return 0;
  }
  
  if (polyData.GetPointer()==NULL)
  {
    // TODO: Error handling
  }
  
  // Points
  igtl::PolyDataPointArray::Pointer pointsArray = polyDataMsg->GetPoints();
  int npoints = pointsArray->GetNumberOfPoints();
  if (npoints > 0)
  {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < npoints; i ++)
    {
      igtlFloat32 point[3];
      pointsArray->GetPoint(i, point);
      points->InsertNextPoint(point); // TODO: use the id returned by this call?
    }
    polyData->SetPoints(points);
  }
  else
  {
    // ERROR: No points defined
  }
  
  // Vertices
  igtl::PolyDataCellArray::Pointer verticesArray =  polyDataMsg->GetVertices();
  int nvertices = verticesArray.IsNotNull() ? verticesArray->GetNumberOfCells() : 0;
  if (nvertices > 0)
  {
    vtkSmartPointer<vtkCellArray> vertCells = vtkSmartPointer<vtkCellArray>::New();
    for (int i = 0; i < nvertices; i ++)
    {
      vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
      
      std::list<igtlUint32> cell;
      verticesArray->GetCell(i, cell);
      //for (unsigned int j = 0; j < cell.size(); j ++) // TODO: is cell.size() always 1?
      //{
      std::list<igtlUint32>::iterator iter;
      iter = cell.begin();
      vertex->GetPointIds()->SetId(i, *iter);
      //}
      vertCells->InsertNextCell(vertex);
    }
    polyData->SetVerts(vertCells);
  }
  
  // Lines
  igtl::PolyDataCellArray::Pointer linesArray = polyDataMsg->GetLines();
  int nlines = linesArray.IsNotNull() ? linesArray->GetNumberOfCells() : 0;
  if (nlines > 0)
  {
    vtkSmartPointer<vtkCellArray> lineCells = vtkSmartPointer<vtkCellArray>::New();
    for(int i = 0; i < nlines; i++)
    {
      vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
      
      std::list<igtlUint32> cell;
      linesArray->GetCell(i, cell);
      polyLine->GetPointIds()->SetNumberOfIds(cell.size());
      std::list<igtlUint32>::iterator iter;
      int j = 0;
      for (iter = cell.begin(); iter != cell.end(); iter ++)
      {
        polyLine->GetPointIds()->SetId(j, *iter);
        j++;
      }
      lineCells->InsertNextCell(polyLine);
    }
    polyData->SetLines(lineCells);
  }
  
  // Polygons
  igtl::PolyDataCellArray::Pointer polygonsArray = polyDataMsg->GetPolygons();
  int npolygons = polygonsArray.IsNotNull() ? polygonsArray->GetNumberOfCells() : 0;
  if (npolygons > 0)
  {
    vtkSmartPointer<vtkCellArray> polygonCells = vtkSmartPointer<vtkCellArray>::New();
    for(int i = 0; i < npolygons; i++)
    {
      vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
      
      std::list<igtlUint32> cell;
      polygonsArray->GetCell(i, cell);
      polygon->GetPointIds()->SetNumberOfIds(cell.size());
      std::list<igtlUint32>::iterator iter;
      int j = 0;
      for (iter = cell.begin(); iter != cell.end(); iter ++)
      {
        polygon->GetPointIds()->SetId(j, *iter);
        j++;
      }
      polygonCells->InsertNextCell(polygon);
    }
    polyData->SetPolys(polygonCells);
  }
  
  // Triangle Strips
  igtl::PolyDataCellArray::Pointer triangleStripsArray = polyDataMsg->GetTriangleStrips();
  int ntstrips = triangleStripsArray.IsNotNull() ? triangleStripsArray->GetNumberOfCells() : 0;
  if (ntstrips > 0)
  {
    vtkSmartPointer<vtkCellArray> tstripCells = vtkSmartPointer<vtkCellArray>::New();
    for(int i = 0; i < ntstrips; i++)
    {
      vtkSmartPointer<vtkTriangleStrip> tstrip = vtkSmartPointer<vtkTriangleStrip>::New();
      
      std::list<igtlUint32> cell;
      triangleStripsArray->GetCell(i, cell);
      tstrip->GetPointIds()->SetNumberOfIds(cell.size());
      std::list<igtlUint32>::iterator iter;
      int j = 0;
      for (iter = cell.begin(); iter != cell.end(); iter ++)
      {
        tstrip->GetPointIds()->SetId(j, *iter);
        j++;
      }
      tstripCells->InsertNextCell(tstrip);
    }
    polyData->SetStrips(tstripCells);
  }
  
  // Attribute
  int nAttributes = polyDataMsg->GetNumberOfAttributes();
  for (int i = 0; i < nAttributes; i ++)
  {
    igtl::PolyDataAttribute::Pointer attribute;
    attribute = polyDataMsg->GetAttribute(i);
    
    vtkSmartPointer<vtkFloatArray> data =
    vtkSmartPointer<vtkFloatArray>::New();
    
    data->SetName(attribute->GetName()); //set the name of the value
    int n = attribute->GetSize();
    
    // NOTE: Data types for POINT (igtl::PolyDataMessage::POINT_*) and CELL
    // (igtl::PolyDataMessage::CELL_*) have the same lower 4 bit.
    // By masking the values with 0x0F, attribute types (either SCALAR, VECTOR, NORMAL,
    // TENSOR, or RGBA) can be obtained. On the other hand, by masking the value
    // with 0xF0, data types (POINT or CELL) can be obtained.
    // See, igtlPolyDataMessage.h in the OpenIGTLink library.
    switch (attribute->GetType() & 0x0F)
    {
      case igtl::PolyDataAttribute::POINT_SCALAR:
      {
        data->SetNumberOfComponents(1);
        break;
      }
      case igtl::PolyDataAttribute::POINT_VECTOR:
      case igtl::PolyDataAttribute::POINT_NORMAL:
      {
        data->SetNumberOfComponents(3);
        break;
      }
      case igtl::PolyDataAttribute::POINT_TENSOR:
      {
        data->SetNumberOfComponents(9); // TODO: Is it valid in Slicer?
        break;
      }
      case igtl::PolyDataAttribute::POINT_RGBA:
      {
        data->SetNumberOfComponents(4); // TODO: Is it valid in Slicer?
        break;
      }
      default:
      {
        // ERROR
        break;
      }
    }
    data->SetNumberOfTuples(n);
    attribute->GetData(static_cast<igtlFloat32*>(data->GetPointer(0)));
    
    if ((attribute->GetType() & 0xF0) == 0) // POINT
    {
      polyData->GetPointData()->AddArray(data);
    }
    else // CELL
    {
      polyData->GetCellData()->AddArray(data);
    }
  }
  // Points RGB
  //std::vector<igtlUint8> pointsRGB = polyDataMsg->GetPointsRGB();
  //int nPointsRGB = polyDataMsg->GetNumberOfPointsRGB();
  if (npoints > 0)
  {
    for (int i = 0; i < npoints; i ++)
    {
      unsigned char color[3] = {255, 0,0};
      colors->InsertNextTupleValue(color);
    }
    polyData->GetPointData()->SetScalars(colors);
  }
  else
  {
    // ERROR: No points defined
  }
  
  polyData->Modified();
  vertexFilter->SetInputData(polyData);
  vertexFilter->Update();
  polyData->ShallowCopy(vertexFilter->GetOutput());
  return true;
  
}


void ConnectionThread()
{
  char*  hostname = "localhost";
  int    port     = 18944;
  
  //------------------------------------------------------------
  // Establish Connection
  
  igtl::ClientSocket::Pointer socket;
  socket = igtl::ClientSocket::New();
  int r = socket->ConnectToServer(hostname, port);
  
  if (r != 0)
  {
    std::cerr << "Cannot connect to the server." << std::endl;
    exit(0);
  }
  
  //------------------------------------------------------------
  // Create a message buffer to receive header
  igtl::MessageHeader::Pointer headerMsg;
  headerMsg = igtl::MessageHeader::New();
  
  //------------------------------------------------------------
  // Allocate a time stamp
  igtl::TimeStamp::Pointer ts;
  ts = igtl::TimeStamp::New();

  while (1)
  {
    // Initialize receive buffer
    headerMsg->InitPack();
    
    // Receive generic header from the socket
    int r = socket->Receive(headerMsg->GetPackPointer(), headerMsg->GetPackSize());
    if (r == 0)
    {
      socket->CloseSocket();
      exit(0);
    }
    if (r != headerMsg->GetPackSize())
    {
      continue;
    }
    
    // Deserialize the header
    headerMsg->Unpack();
    
    // Get time stamp
    igtlUint32 sec;
    igtlUint32 nanosec;
    
    headerMsg->GetTimeStamp(ts);
    ts->GetTimeStamp(&sec, &nanosec);
    
    
    // Check data type and receive data body
    if (strcmp(headerMsg->GetDeviceType(), "POLYDATA") == 0)
    {
      bool success = ReceivePolyDataStream(socket, headerMsg);
      if (success)
      {
        mapper->SetInputData(polyData);
        renderer->GetRenderWindow()->Render();
      }
    }
    else
    {
      std::cerr << "Receiving : " << headerMsg->GetDeviceType() << std::endl;
      socket->Skip(headerMsg->GetBodySizeToRead(), 0);
    }
  }
  
  //------------------------------------------------------------
  // Close connection (The example code never reaches this section ...)
  
  socket->CloseSocket();
}

class customMouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
  static customMouseInteractorStyle* New();
  vtkTypeMacro(customMouseInteractorStyle, vtkInteractorStyleTrackballCamera);
  virtual void OnLeftButtonDown()
  {
    std::cout << "Pressed left mouse button." << std::endl;
    // Forward events
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    interactionActive = true;
  }
  
  virtual void OnLeftButtonUp()
  {
    std::cout << "Release left mouse button." << std::endl;
    // Forward events
    vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
    interactionActive = false;
    conditionVar->Signal();
  }
  
  virtual void OnMiddleButtonDown()
  {
    std::cout << "Pressed middle mouse button." << std::endl;
    // Forward events
    vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
  }
  
  virtual void OnRightButtonDown()
  {
    std::cout << "Pressed right mouse button." << std::endl;
    // Forward events
    vtkInteractorStyleTrackballCamera::OnRightButtonDown();
  }
};
vtkStandardNewMacro(customMouseInteractorStyle);

int main(int argc, char* argv[])
{
  conditionVar = igtl::ConditionVariable::New();
  localMutex = igtl::SimpleMutexLock::New();
  interactionActive = false;
  // Create a mapper
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(polyData);
  mapper->SetColorModeToDefault();
  
  // Create an actor
  actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  
  // A renderer and render window
  
  // Create a renderer, render window, and interactor
  renderer = vtkSmartPointer<vtkRenderer>::New();
  renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  
  // An interactor
  renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  
  // Add the actors to the scene
  renderer->AddActor(actor);
  renderer->SetBackground(0.2, 0.3, 0.4);
  
  // Render an image (lights and cameras are created automatically)
  
  renderWindow->Render();
  igtl::MultiThreader::Pointer threadViewer = igtl::MultiThreader::New();
  // Begin mouse interaction
  colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetName("Color");
  colors->SetNumberOfComponents(3);
  vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
  cloud = vtkSmartPointer<vtkPoints>::New();
  polyData = vtkSmartPointer<vtkPolyData>::New();
  
  mapper->SetInputData(polyData);
  // Render an image (lights and cameras are created automatically)
  renderWindow->SetSize(1000, 600);
  renderWindow->Render();
  vtkInteractorObserver *currentStyle =
  renderWindowInteractor->GetInteractorStyle();
  std::cout << "currentStyle class name: " << currentStyle->GetClassName() << std::endl;
  
  vtkInteractorStyleSwitch *iss = vtkInteractorStyleSwitch::SafeDownCast(currentStyle);
  vtkInteractorObserver *actualStyle = iss->GetCurrentStyle();
  std::cout << "actualStyle class name: " <<
  actualStyle->GetClassName() << std::endl;
  vtkSmartPointer<customMouseInteractorStyle> style = vtkSmartPointer<customMouseInteractorStyle>::New();
  renderWindowInteractor->SetInteractorStyle( style );
  std::cout << "currentStyle class name: " << renderWindowInteractor->GetInteractorStyle()->GetClassName() << std::endl;
  threadViewer->SpawnThread((igtl::ThreadFunctionType) &ConnectionThread, NULL);
  renderWindowInteractor->Start();

  while(1)
  {
    igtl::Sleep(1000);
  }
}

