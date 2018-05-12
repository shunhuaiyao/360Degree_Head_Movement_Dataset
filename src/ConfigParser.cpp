// Author: Xavier Corbillon
// IMT Atlantique
#include "ConfigParser.hpp"

#include "MeshCube.hpp"
#include "MeshEquiAngularCube.hpp"
#include "MeshCubeEquiUV.hpp"
#include "MeshEqualarea.hpp"
#include "MeshSegmentedSphere.hpp"
#include "ShaderTextureStatic.hpp"
#include "ShaderTextureVideo.hpp"
#include "LogWriter.hpp"
#include "PublisherLogMQ.hpp"

//Library includes
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

//standard library
#include <iostream>
#include <stdexcept>

using namespace IMT;


void ConfigParser::Init(void)
{
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(m_pathToConfigFile, pt);

  std::cout << "Start parsing " << m_pathToConfigFile << " configuration file\n";
  //Read the Config section
  auto textureConfig = pt.get<std::string>("Config.textureConfig");
  auto projectionConfig = pt.get<std::string>("Config.projectionConfig");
  auto logWriterConfig = pt.get<std::string>("Config.logWriterConfig");
  auto publisherLogConfig = pt.get<std::string>("Config.publisherLogConfig");

  std::cout << "Parse the texture configuration: section "<< textureConfig <<"\n";

  auto textureType = pt.get<std::string>(textureConfig+".type");
  //Static picture
  if (textureType == "static")
  {
    auto pathToPicture = pt.get<std::string>(textureConfig+".pathToPicture");
    m_outputShaderTexture = std::make_shared<ShaderTextureStatic>(pathToPicture);
  }
  //Video
  else if (textureType == "video")
  {
    auto pathToVideo = pt.get<std::string>(textureConfig+".pathToVideo");
    size_t nbFrame = pt.get<size_t>(textureConfig+".nbFrame");
    size_t  bufferSize = pt.get<size_t>(textureConfig+".bufferSize");
    float startOffsetInSecond = pt.get<float>(textureConfig+".startOffsetInSecond");
    std::cout << "bufferSize " <<bufferSize << std::endl;
    m_outputShaderTexture = std::make_shared<ShaderTextureVideo>(pathToVideo, nbFrame, bufferSize, startOffsetInSecond);
  }
  else
  {
    throw(std::invalid_argument("Not supported texture type: "+textureType));
  }

  std::cout << "Parse the projection configuration: section "<< projectionConfig <<"\n";

  auto projectionType = pt.get<std::string>(projectionConfig+".type");
  if (projectionType == "CubeMap" || projectionType == "AdjustedCubeMap")
  {
    m_outputMesh = std::make_shared<MeshCube>(10.0f, 6*2*30*30);
  }
  else if (projectionType == "EquiAngularCubeMap")
  {
	m_outputMesh = std::make_shared<MeshEquiAngularCube>(10.0f, 6*2*30*30);
  }
  else if (projectionType == "Equirectangular")
  {
    m_outputMesh = std::make_shared<MeshCubeEquiUV>(10.0f, 6*2*30*30);
  }
  else if (projectionType == "Equalarea")
  {
    m_outputMesh = std::make_shared<MeshEqualarea>(10.0f, 6*2*30*30);
  }
  else if (projectionType == "SegmentedSphere")
  {
    m_outputMesh = std::make_shared<MeshSegmentedSphere>(10.0f, 6*2*1040*1040);
  }
  else
  {
    throw(std::invalid_argument("Not supported projection type: "+projectionType));
  }

  std::cout << "Parse the log writer configuration: section "<< logWriterConfig <<"\n";

  auto logWriterOutputDirPath = pt.get<std::string>(logWriterConfig+".outputDirPath");
  auto logWriterOutputId = pt.get<std::string>(logWriterConfig+".outputId");
  m_outputLogWriter = std::make_shared<LogWriter>(logWriterOutputDirPath, logWriterOutputId);

  std::cout << "Parse the publisher log configuration: section "<< publisherLogConfig << "\n";

  auto publisherLogPort = pt.get<size_t>(publisherLogConfig+".port");
  m_outputPublisherLogMQ = std::make_shared<PublisherLogMQ>();
  m_outputPublisherLogMQ->Init(publisherLogPort);
}
