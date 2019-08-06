/*
Copyright (C) 2012 Jeffrey M Brown

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with This program. If not, see <http://www.gnu.org/licenses/>.
*/


// config.cc - parses and stores values from configuration files

#include "config.h"

Config::Config() {}

Config::Config(const std::string& filename) {
  load(filename);
}

Config::~Config() {}

void Config::load(const std::string& filename) {
  // open file
  std::ifstream ifs(filename.c_str());
  if(!ifs) throw ConfigFileError("Cannot open input file: " + filename);
  
  // read lines from file
  std::string line;
  std::vector<std::string> text;
  while(!ifs.eof()) {
    getline(ifs, line);
    // ignore comment or blank lines
    if(line.length() > 0 && line[0] != '#') 
      text.push_back(trimWhitespace(line));
  }

  // iterate through text and build the dictionary
  std::string sectionName, key, value, fullPathKey;
  std::vector<std::string>::const_iterator text_iter;
  for(text_iter=text.begin(); text_iter!=text.end(); text_iter++) {
    line = *text_iter; // current line of text
    // check if line is a section header
    if(*line.begin() == '[' && *line.rbegin() == ']') {
      sectionName = trimWhitespace(line.substr(1, line.length()-2));
    }
    else { // line might contain: key = value
      int i = (int)  line.find("=");
      if(i == -1) // failed to find '='
	throw ConfigFileError("Under section '"+sectionName+"', malformed expression line: '"+line+"'\n");

      // split expression by '=' into key and value
      key = trimWhitespace(line.substr(0, i));
      value = trimWhitespace(line.substr(i+1));
      if( (int) !key.length() > 0 || (int) !value.length() > 0) // if key or value contains only whitespace
	throw ConfigFileError("Under section '"+sectionName+"', missing key or value: '"+key+"="+value+"'\n");

      // store keys as 'section/key'
      fullPathKey = sectionName + "/" + key;
      settings[fullPathKey] = value;
    }
  }
}

// save settings map back to a configuration file format
void Config::save(const std::string& filename) {
  // open file
  std::ofstream ofs(filename.c_str());
  if(!ofs) throw ConfigFileError("Cannot open output file: " + filename);
  
  // iterate through the dictionary
  std::string fullPathKey, sectionName, key, value, currentSectionName;
  std::map<std::string, std::string>::const_iterator map_iter;
  for(map_iter=settings.begin(); map_iter!=settings.end(); map_iter++) {
    fullPathKey = map_iter->first;
    value = map_iter->second;
    
    // split by '/'
    int i = (int)  fullPathKey.find("/");
    sectionName = fullPathKey.substr(0, i);
    key = fullPathKey.substr(i+1);

    // check if we need to start a new section
    if(currentSectionName == sectionName) { // same section still
      ofs << key << " = " << value << "\n";
    }
    else { // create a new section
      ofs << "\n" << "[" << sectionName << "]\n";
      ofs << key << " = " << value << "\n";
      currentSectionName = sectionName; // update the section we are under
    }
  }
}

// print out all of the key/value pairs of settings
void Config::printSettings(std::ostream& os) {
  os << "*** Settings ***\n";
  std::map<std::string, std::string>::const_iterator map_iter;
  for(map_iter=settings.begin(); map_iter!=settings.end(); map_iter++) {
    os << map_iter->first << ": " << map_iter->second << "\n";
  }
}

// trim leading and trailing whitespace
std::string Config::trimWhitespace(const std::string& str) {
  // find first and last characters that are not spaces
  int first = (int) str.find_first_not_of(" "); // returns -1 if only spaces present
  int last = (int)  str.find_last_not_of(" ");
  
  // check if we found anything other than spaces
  if(first >= 0 && last >= 0) {
    return str.substr(first, last-first+1);
  }
  else { // string only contains spaces
    return std::string("");
  }
}

