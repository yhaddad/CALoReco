#idndef _HYUTILS_H
#define _HYUTILS_H

#include <map>
#include <algorithm>
#include <vector>
#include <marlin/tinyxml.h>

using namespace std;



class HyUtils
{
 public:
  void xml_geom_reader(std::string xmlfile){
    TiXmlDocument doc(xmlfile.c_str());
    bool load_key = doc.LoadFile();  
    if(load_key){
      std::cout << green << "File : " << xmlfile.c_str() << normal <<std::endl;
      // tout ici 
      TiXmlHandle hDoc(&doc);
      TiXmlElement* pElem;
      TiXmlHandle hRoot(0);
      // name block
      {
	pElem=hDoc.FirstChildElement().Element();
	// should always have a valid root but handle gracefully if it does
	if (!pElem) std::cout << red << "error elem" << normal << std::endl;
	std::cout << green << pElem->Value() << normal << std::endl;
	
	// save this for later
	hRoot=TiXmlHandle(pElem);
      }
      // parameters block
      {
	m_parameters.clear();
	pElem=hRoot.FirstChild("parameter").Element();
	std::string key = pElem->Attribute("name");
	std::cout << green << key.c_str() << normal << std::endl; 
#ifdef DEBUG
      {
      	std::cout << green
      		  <<"parameter : " 
      		  << pElem->Attribute("name") 
      		  << normal 
      		  << std::endl;
      }
#endif
      std::vector<std::string> lines;
      {
	stringbuf *pbuf;
	std::string value = pElem->GetText() ;
	std::vector<std::string> lines;
	istringstream iss(value);
	copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter<vector<string> >(lines));
	for(int iline = 0; iline < lines.size(); iline++){
	  std::string line = lines.at(iline);
	  std::cout << red << line << normal << std::endl;
	  
	  stringstream ss( line.c_str() );
	  vector<string> result;
	  
	  LayerID mapp;
	  int Dif_id;
	  while( ss.good() )
	    {
	      string substr;
	      getline( ss, substr, ',' );
	      result.push_back( substr );
	    }
	  istringstream ( result.at(0) ) >> Dif_id;
	  istringstream ( result.at(1) ) >> mapp.K;
	  istringstream ( result.at(2) ) >> mapp.DifX;
	  istringstream ( result.at(3) ) >> mapp.DifY;
	  istringstream ( result.at(4) ) >> mapp.IncX;
	  istringstream ( result.at(5) ) >> mapp.IncY;
	  _mapping[Dif_id] = mapp;
	}
      }
      pElem = pElem->NextSiblingElement();
      // ChamberGeom  Node.
      {
#ifdef DEBUG
	{
	  std::cout << green
		    <<"parameter : " 
		    << pElem->Attribute("name") 
		    << normal 
		    << std::endl;
	}
#endif
	std::vector<std::string> lines;
	{
	  stringbuf *pbuf;
	  std::string value = pElem->GetText() ;
	  std::vector<std::string> lines;
	  istringstream iss(value);
	  copy(istream_iterator<string>(iss),
	       istream_iterator<string>(),
	       back_inserter<vector<string> >(lines));
	  for(int iline = 0; iline < lines.size(); iline++){
	    std::string line = lines.at(iline);
	    std::cout << red << line << normal << std::endl;
	    
	    stringstream ss( line.c_str() );
	    vector<string> result;
	    
	    double position;
	    int Dif_id;
	    while( ss.good() )
	      {
		string substr;
		getline( ss, substr, ',' );
		result.push_back( substr );
	      }
	    istringstream ( result.at(0) ) >> Dif_id;
	    istringstream ( result.at(3) ) >> position;
	    
	    _chamber_pos[Dif_id] = position;
	  }
	}
      }
    }
  }else{
    std::cout << red << "Failed to load file : " << xmlfile.c_str() << normal <<std::endl;
  }
}

};
#endif
