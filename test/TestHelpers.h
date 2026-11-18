#ifndef TestHelpers_h
#define TestHelpers_h

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <atomic>
#include <unistd.h>
#include <functional>
#include <csignal>
#include <fcntl.h>
#include <sys/wait.h>
#include "../src/AlleleSpecificExpression.h"

namespace TestHelpers
{

inline CellMatch makeMatch(const std::string& cell, const std::string& variant, bool alt, size_t count)
{
	CellMatch result;
	result.cell = cell;
	result.variant = variant;
	result.alt = alt;
	result.count = count;
	return result;
}

// A temporary file which is removed when the object goes out of scope.
class TempFile
{
public:
	explicit TempFile(const std::string& suffix = "")
	{
		static std::atomic<size_t> counter { 0 };
		path = (std::filesystem::temp_directory_path() / ("xcieem_unittest_" + std::to_string(getpid()) + "_" + std::to_string(counter++) + suffix)).string();
	}
	TempFile(const std::string& suffix, const std::string& contents) : TempFile(suffix)
	{
		write(contents);
	}
	~TempFile()
	{
		std::error_code ec;
		std::filesystem::remove(path, ec);
	}
	TempFile(const TempFile&) = delete;
	TempFile& operator=(const TempFile&) = delete;
	void write(const std::string& contents) const
	{
		std::ofstream file { path, std::ios::binary };
		file << contents;
	}
	std::string read() const
	{
		std::ifstream file { path, std::ios::binary };
		std::stringstream sstr;
		sstr << file.rdbuf();
		return sstr.str();
	}
	std::vector<std::string> readLines() const
	{
		std::vector<std::string> result;
		std::ifstream file { path };
		std::string line;
		while (std::getline(file, line)) result.push_back(line);
		return result;
	}
	std::string path;
};

// Runs the function in a forked child process and returns true if the child was killed by SIGABRT.
// Used to test the error paths which call std::abort(). The child's stderr is discarded.
inline bool abortsInChild(const std::function<void()>& function)
{
	fflush(stdout);
	fflush(stderr);
	pid_t pid = fork();
	if (pid == 0)
	{
		// doctest installs its own crash handlers; the child should just die on abort
		signal(SIGABRT, SIG_DFL);
		int devnull = open("/dev/null", O_WRONLY);
		if (devnull >= 0)
		{
			dup2(devnull, 1);
			dup2(devnull, 2);
		}
		function();
		_exit(0);
	}
	int status = 0;
	waitpid(pid, &status, 0);
	return WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
}

inline std::vector<std::string> lines(const std::string& text)
{
	std::vector<std::string> result;
	std::istringstream stream { text };
	std::string line;
	while (std::getline(stream, line)) result.push_back(line);
	return result;
}

}

#endif
