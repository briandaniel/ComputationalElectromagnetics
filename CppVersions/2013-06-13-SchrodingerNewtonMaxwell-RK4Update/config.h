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


// config.h - parses and stores values from configuration files

/*
Configuration files should have the following format:
--------------------
[Pump]
type = Gaussian
width = 20

[Atom]
B = 1.05
--------------------
where key/value pairs can be grouped under sections, denoted
by square braces.

These pairs will be stored in a map<string, string> called settings
with each section prepended to the individual keys listed below it.

Pump/type: Gaussian
Pump/width: 20
Atom/B: 1.05

Each value can be read out as whatever type is passed 
into the get method thanks to the magic of templates.

Config cf(configfile);
double B;
cf.get("Atom/B", B) // B is now 1.05

or

int B;
cf.get("Atom/B", B) // B is now 1

*/

#ifndef CONFIG_H_
#define CONFIG_H_

#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <iostream>

// general exception
struct ConfigFileError : std::runtime_error {
  ConfigFileError(const std::string& message)
    :std::runtime_error(message) {}
};

// exception to be thrown when a map is missing a key
struct KeyError : std::runtime_error {
  KeyError(const std::string& key)
    :std::runtime_error("Could not find key: '" + key + "'") {}
};


// class to parse configuration files and store their contents in a map
class Config {
public:
  Config();
  Config(const std::string& filename);    // loads file upon instantiation
  virtual ~Config();

  void load(const std::string& filename); // add this file to settings
  void save(const std::string& filename); // save the contents of settings to a file
  
  // main method of getting the value of a key from settings
  template <class T>
  void get(const std::string& key, T& value);

  // main method for setting the value of a key from settings
  template <class T>
  void set(const std::string& key, const T& value);
  
  // print out the keys and values
  void printSettings(std::ostream& os);

private:
  std::map<std::string, std::string> settings;
  std::stringstream ss; // used for type conversions
  
  // removes leading and trailing whitespace
  std::string trimWhitespace(const std::string& str);
};


// templates must be defined in header to avoid linking error
template <class T>
void Config::get(const std::string& key, T& value) {
  std::map<std::string, std::string>::iterator iter = settings.find(key);
  if(iter != settings.end()) { // key exists
    ss << iter->second;    // send string value into stream
    ss >> value;           // read out value into proper type
    
    ss.str(std::string()); // set contents to an empty string
    ss.clear();            // clear any errors (most likely eof)
  }
  else throw KeyError(key);
}

template <class T>
void Config::set(const std::string& key, const T& value) {
  ss << value;              // read in value
  settings[key] = ss.str(); // save string value for key

  ss.str(std::string());    // set contents to an empty string
  ss.clear();               // clear any errors (most likely eof)
}



#endif // CONFIG_H_
