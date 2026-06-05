#ifndef Common_h
#define Common_h

#include <string>
#include <vector>
#include <tuple>
#include <unordered_map>

std::string matHapName(const bool phasesAreMatPat);
std::string patHapName(const bool phasesAreMatPat);
std::vector<std::string> split(const std::string& raw, const char separator);
std::string join(const char separator, const std::vector<std::string>& strings);
std::string lowercase(std::string raw);
std::vector<std::string> getVariantOrder(const std::unordered_map<std::string, size_t>& nameToIndex);
std::vector<std::string> getCellOrder(const std::unordered_map<std::string, size_t>& nameToIndex);
size_t parseVariantPosition(const std::string& name);
std::tuple<std::string, size_t, size_t> parseBedRegion(const std::string& region);
double getBinomialPValueGreaterThan(const double p, const size_t successes, const size_t trials);
double getBinomialPValueLessThan(const double p, const size_t successes, const size_t trials);
bool hasExtension(const std::string& filename, const std::string& extension);

#endif
