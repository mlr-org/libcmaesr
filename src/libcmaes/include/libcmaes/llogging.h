/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Copyright (c) 2014 Inria
 * Author: Emmanuel Benazera <emmanuel.benazera@lri.fr>
 *
 * This file is part of libcmaes.
 *
 * libcmaes is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libcmaes is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with libcmaes.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LLOGGING_H
#define LLOGGING_H

#include <libcmaes/libcmaes_config.h>

#if defined(HAVE_GLOG)

#include <glog/logging.h>

#elif defined(LIBCMAES_LOG_TO_R)

#include <ostream>
#include <streambuf>
#include <string>

extern "C" void Rprintf(const char *format, ...);

namespace libcmaes {
static std::string INFO = "INFO";
static std::string WARNING = "WARNING";
static std::string ERROR = "ERROR";
static std::string FATAL = "FATAL";

static std::ostream nullstream(0);

// Minimal streambuf that forwards output directly to Rprintf (unbuffered)
class RPrintfBuf : public std::streambuf {
public:
  using int_type = std::streambuf::int_type;
  using traits_type = std::streambuf::traits_type;

protected:
  // Write a block of characters. iostreams prefer this path for operator<< of strings
  // and formatted output. We forward the entire chunk to Rprintf in one call.
  // Returning 'n' tells iostreams the full block was consumed.
  std::streamsize xsputn(const char *s, std::streamsize n) override {
    if (n > 0) Rprintf("%.*s", static_cast<int>(n), s);
    return n;
  }

  // Write a single character (fallback used by iostreams when xsputn isn't used).
  // Called for individual chars, including from operator<< for char and as a
  // fallback in some formatted paths. We emit that character via Rprintf.
  // Must return the written char on success, or traits::not_eof on EOF.
  int_type overflow(int_type ch) override {
    // Only process and emit when ch is not the EOF marker (overflow may be called with EOF just to flush)
    if (!traits_type::eq_int_type(ch, traits_type::eof())) {
      const char c = traits_type::to_char_type(ch);
      Rprintf("%c", c);
      return ch;
    }
    return traits_type::not_eof(ch);
  }

  // Flush any buffered content. This minimal implementation is unbuffered,
  // so there’s nothing to flush; return 0 to signal success so std::endl works.
  int sync() override { return 0; }
};

// Global ostream that uses RPrintfBuf
static RPrintfBuf rprintf_buf;
static std::ostream rcout(&rprintf_buf);

inline std::ostream &LOG(const std::string &severity, std::ostream &out = rcout) {
  out << severity << " - ";
  return out;
}

inline std::ostream &LOG_IF(const std::string &severity, const bool &condition, std::ostream &out = rcout) {
  if (condition)
    return LOG(severity, out);
  else
    return nullstream;
}
} // namespace libcmaes

#else

#include <iostream>

namespace libcmaes {
static std::string INFO = "INFO";
static std::string WARNING = "WARNING";
static std::string ERROR = "ERROR";
static std::string FATAL = "FATAL";

static std::ostream nullstream(0);

inline std::ostream &LOG(const std::string &severity, std::ostream &out = std::cout) {
  out << severity << " - ";
  return out;
}

inline std::ostream &LOG_IF(const std::string &severity, const bool &condition, std::ostream &out = std::cout) {
  if (condition)
    return LOG(severity, out);
  else
    return nullstream;
}
} // namespace libcmaes

#endif

#endif
