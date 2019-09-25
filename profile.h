#pragma once

class LogDuration {
 public:
  explicit LogDuration(const std::string& msg = "")
      : message(msg + ": "), start(std::chrono::steady_clock::now()) {
  }

  int getTime() {
    auto finish = std::chrono::steady_clock::now();
    auto dur = finish - start;
    return std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
  }

  ~LogDuration() {
    auto finish = std::chrono::steady_clock::now();
    auto dur = finish - start;
    std::ostringstream os;
    os << message
       << std::chrono::duration_cast<std::chrono::milliseconds>(dur).count()
       << " ms" << std::endl;
    std::cerr << os.str() << std::endl;
  }
 private:
  std::string message;
  std::chrono::steady_clock::time_point start;
};

#ifndef UNIQ_ID
#define UNIQ_ID_IMPL(lineno) _a_local_var_##lineno
#define UNIQ_ID(lineno) UNIQ_ID_IMPL(lineno)
#endif

#define LOG_DURATION(message) \
  LogDuration UNIQ_ID(__LINE__){message};
