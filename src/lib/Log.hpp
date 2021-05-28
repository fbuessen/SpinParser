/**
 * @file Log.hpp
 * @author Finn Lasse Buessen
 * @brief Lightweight logging interface with output filtering.
 *
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <iostream>
#include <functional>
#include <boost/date_time.hpp>


namespace Log
{
	/**
	 * @brief Specify log level for output filtering. Any Log level includes also output from lower output levels.
	 */
	enum struct LogLevel
	{
		None = 0, ///< No output. 
		Error = 1, ///< Print error messages. 
		Warning = 2, ///< Print warnings. 
		Info = 3, ///< Print information. 
		Debug = 4 ///< Print debug output. 
	};

	class Logstream;

	/**
	 * @brief Logstream manipulator type with internal structure. 
	 */
	typedef std::function<Logstream & (Logstream &)> StructManipulator;

	/**
	 * @brief Logstream manipulator type without internal structure.
	 */
	typedef Logstream &(*Manipulator)(Logstream &);

	/**
	 * @brief Log stream object for simple output filtering.
	 * @details The Logstream object provides output filtering accordign to a selected log level.
	 * The filtered output is written to stdout.
	 * Any output generated is printed with a timestamp, measring the time since creation of the Logstream object.
	 */
	class Logstream
	{
		friend Logstream &endl(Logstream &ls);
		friend StructManipulator setDisplayLogLevel(const Log::LogLevel logLevel);
		friend StructManipulator setLogLevel(const Log::LogLevel logLevel);
	public:
		/**
		 * @brief Construct a new Logstream object with log filtering to Log::LogLevel::Info.
		 */
		Logstream(std::ostream &logTarget) : _logTarget(logTarget), _constructionTime(boost::posix_time::microsec_clock::local_time())
		{
			_carriageReturn = true;
			_displayLogLevel = Log::LogLevel::Info;
			_streamLogLevel = Log::LogLevel::Info;
		};

		/**
		 * @brief Output operator for messages of arbitrary type. Will accept any log object that implements the output operator for stdout.
		 *
		 * @tparam Type of the log message.
		 * @param ls Target Logstream object.
		 * @param rhs Log message to be written.
		 * @return Logstream& Reference to Logstream operand.
		 */
		template <class T> friend Logstream &operator<<(Logstream &ls, const T &rhs)
		{
			if (ls._displayLogLevel < ls._streamLogLevel) return ls;
			if (ls._streamLogLevel == Log::LogLevel::Warning) ls._logTarget << "\033[31m";
			if (ls._carriageReturn)
			{
				float dt = float((boost::posix_time::microsec_clock::local_time() - ls._constructionTime).total_microseconds());
				ls._logTarget << "[" << std::fixed << std::setprecision(6) << dt / 1000000.0f << "][";
				if (ls._streamLogLevel == LogLevel::Debug) ls._logTarget << "D";
				else if (ls._streamLogLevel == LogLevel::Info) ls._logTarget << "I";
				else if (ls._streamLogLevel == LogLevel::Warning) ls._logTarget << "W";
				else if (ls._streamLogLevel == LogLevel::Error) ls._logTarget << "E";
				ls._logTarget << "] ";
				ls._carriageReturn = false;
			}
			ls._logTarget << rhs;
			if (ls._streamLogLevel == Log::LogLevel::Warning) ls._logTarget << "\033[0m";

			return ls;
		}

		/**
		 * @brief Output operator for LogLevel objects to set filter level.
		 *
		 * @param rhs New LogLevel value.
		 * @return Logstream& Reference to self.
		 */
		Logstream &operator<<(const Log::LogLevel &rhs)
		{
			_streamLogLevel = rhs;
			return *this;
		}

		/**
		 * @brief Apply plain output manipulator to log stream.
		 *
		 * @param rhs Logstream manipulator.
		 * @return Logstream& Reference to self.
		 */
		Logstream &operator<<(const Manipulator &rhs)
		{
			return rhs(*this);
		}

		/**
		 * @brief Apply structured output manipulator to log stream.
		 *
		 * @param rhs Logstream manipulator.
		 * @return Logstream& Reference to self.
		 */
		Logstream &operator<<(const StructManipulator &rhs)
		{
			return rhs(*this);
		}

	private:
		std::ostream &_logTarget; ///< Output stream. 
		boost::posix_time::ptime _constructionTime; ///< Creation time of the Logstream object. 

		bool _carriageReturn; ///< If set to true, the next log message is printed with timestamp. 
		Log::LogLevel _displayLogLevel; ///< Selected output filtering. @see Log::LogLevel
		Log::LogLevel _streamLogLevel; ///< Current log level for new log messages. 
	};

	#pragma region Manipulators
	/**
	 * @brief Output modifier to print new line.
	 *
	 * @param ls Target Logstream object.
	 * @return Logstream& Reference to self.
	 */
	inline Logstream &endl(Logstream &ls)
	{
		if (ls._displayLogLevel >= ls._streamLogLevel)
		{
			ls._logTarget << std::endl;
			ls._carriageReturn = true;
		}
		return ls;
	}

	/**
	 * @brief Generate output modifier to set the display log level filter.
	 *
	 * @param logLevel Filter level.
	 * @return StructManipulator Logstream manipulator.
	 */
	inline StructManipulator setDisplayLogLevel(const Log::LogLevel logLevel)
	{
		auto manipulator = [logLevel](Logstream &ls)->Logstream &
		{
			ls._displayLogLevel = logLevel;
			return ls;
		};
		return manipulator;
	}

	/**
	 * @brief Generate output modifier to change the log level.
	 *
	 * @param logLevel New log level.
	 * @return StructManipulator Logstream manipulator.
	 */
	inline StructManipulator setLogLevel(const Log::LogLevel logLevel)
	{
		auto manipulator = [logLevel](Logstream &ls)->Logstream &
		{
			ls._streamLogLevel = logLevel;
			return ls;
		};
		return manipulator;
	}
	#pragma endregion

	extern Logstream log; ///< Global Logstream instance. 
}