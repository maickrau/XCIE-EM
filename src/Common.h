#ifndef Common_h
#define Common_h

#include <string>
#include <vector>
#include <tuple>
#include <unordered_map>

std::string matHapName(const bool phasesAreMatPat);
std::string patHapName(const bool phasesAreMatPat);
std::vector<std::string> split(const std::string& raw, const char separator);
std::string lowercase(std::string raw);
std::vector<std::string> getVariantOrder(const std::unordered_map<std::string, size_t>& nameToIndex);
std::vector<std::string> getCellOrder(const std::unordered_map<std::string, size_t>& nameToIndex);
size_t parseVariantPosition(const std::string& name);
std::tuple<std::string, size_t, size_t> parseBedRegion(const std::string& region);

#endif
