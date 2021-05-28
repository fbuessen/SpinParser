/**
 * @file Integrator.hpp
 * @author Finn Lasse Buessen
 * @brief One-dimensional numerical integration routines. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <functional>
#include "lib/Assert.hpp"
#include "FrgCommon.hpp"

namespace Integrator
{
	/**
	 * @brief One-dimensional trapezoidal integration in one-dimensional frequency space. 
	 * 
	 * @tparam T Integrand type. 
	 * @param min Iterator to lower boundary frequency value. 
	 * @param max Iterator to upper boundary frequency value. 
	 * @param integrand Integrand function. 
	 * @return T Value of the integral. 
	 */
	template <class T> T integrate(const FrequencyIterator min, const FrequencyIterator max, const std::function<T(float)> &integrand)
	{
		ASSERT(*min <= *max, "Lower integration boundary must not be larger than upper boundary. ");

		FrequencyIterator umin(min);
		if (umin != max)
		{
			T integral = (*(umin + 1) - *umin) * integrand(*umin);
			while (++umin != max) integral += (*(umin + 1) - *(umin - 1)) * integrand(*umin);
			integral += (*umin - *(umin - 1)) * integrand(*umin);
			return 0.5f * integral;
		}
		else return 0.0f;
	}

	/**
	 * @brief One-dimensional trapezoidal integration in one-dimensional frequency space. 
	 * 
	 * @tparam T Integrand type. 
	 * @param min Lower boundary frequency value. 
	 * @param max Iterator to upper boundary frequency value. 
	 * @param integrand Integrand function. 
	 * @return T Value of the integral. 
	 */
	template <class T> T integrateWithObscureLeftBoundary(const float min, const FrequencyIterator max, const std::function<T(float)> &integrand)
	{
		ASSERT(min <= *max, "Lower integration boundary must not be larger than upper boundary. ");

		FrequencyIterator umin = FrgCommon::frequency().greater(min);
		if (umin != max)
		{
			T integral = (*umin - min) * integrand(min);
			integral += (*(umin + 1) - min) * integrand(*umin);
			while (++umin != max) integral += (*(umin + 1) - *(umin - 1)) * integrand(*umin);
			integral += (*umin - *(umin - 1)) * integrand(*umin);
			return 0.5f * integral;
		}
		else return 0.5f * (*umin - min) * (integrand(min) + integrand(*umin));
	}

	/**
	 * @brief One-dimensional trapezoidal integration in one-dimensional frequency space. 
	 * 
	 * @tparam T Integrand type. 
	 * @param min Iterator to lower boundary frequency value. 
	 * @param max Upper boundary frequency value. 
	 * @param integrand Integrand function. 
	 * @return T Value of the integral. 
	 */
	template <class T> T integrateWithObscureRightBoundary(const FrequencyIterator min, const float max, const std::function<T(float)> &integrand)
	{
		ASSERT(*min <= max, "Lower integration boundary must not be larger than upper boundary. ");

		FrequencyIterator umin(min);
		FrequencyIterator umax = FrgCommon::frequency().lesser(max);
		if (umin != umax)
		{
			T integral = (*(umin + 1) - *umin) * integrand(*umin);
			while (++umin != umax) integral += (*(umin + 1) - *(umin - 1)) * integrand(*umin);
			integral += (max - *(umin - 1)) * integrand(*umin);
			integral += (max - *umin) * integrand(max);
			return 0.5f * integral;
		}
		else return 0.5f * (max - *umin) * (integrand(max) + integrand(*umin));
	}

	/**
	 * @brief One-dimensional trapezoidal integration in one-dimensional frequency space. 
	 * 
	 * @tparam T Integrand type. 
	 * @param min Lower boundary frequency value. 
	 * @param max Upper boundary frequency value. 
	 * @param integrand Integrand function. 
	 * @return T Value of the integral. 
	 */
	template <class T>
	T integrateWithObscureBoundaries(const float min, const float max, const std::function<T(float)> &integrand)
	{
		ASSERT(min <= max, "Lower integration boundary must not be larger than upper boundary. ");

		FrequencyIterator umin = FrgCommon::frequency().greater(min);
		FrequencyIterator umax = FrgCommon::frequency().lesser(max);
		if (umax >= umin)
		{
			T integral = (*umin - min) * integrand(min);
			if (umax != umin)
			{
				integral += (*(umin + 1) - min) * integrand(*umin);
				while (++umin != umax) integral += (*(umin + 1) - *(umin - 1)) * integrand(*umin);
				integral += (max - *(umin - 1)) * integrand(*umin);
			}
			else integral += (max - min) * integrand(*umin);
			integral += (max - *umin) * integrand(max);
			return 0.5f * integral;
		}
		else return 0.5f * (max - min) * (integrand(max) + integrand(min));
	}
}

namespace ImplicitIntegrator
{
	/**
	 * @brief One-dimensional trapezoidal integration in one-dimensional frequency space. 
	 * First argument of the integrand function is the frequency value, the second argument is the return value of the integrand, passed by reference. 
	 * 
	 * @tparam T Integrand type. 
	 * @param[in] min Lower boundary frequency value. 
	 * @param[in] max Iterator to upper boundary frequency value. 
	 * @param[in] integrand Integrand function. 
	 * @param[out] integrandBuffer Return value buffer of the integrand. 
	 * @param[out] resultBuffer Value of the integral. 
	 */
	template <class T>
	void integrateWithObscureLeftBoundary(const float min, const FrequencyIterator max, const std::function<void(float, T &)> &integrand, T &integrandBuffer, T &resultBuffer)
	{
		ASSERT(&integrandBuffer != &resultBuffer);
		ASSERT(min <= *max, "Lower integration boundary must not be larger than upper boundary. ");

		resultBuffer.reset();
		FrequencyIterator umin = FrgCommon::frequency().greater(min);
		if (umin != max)
		{
			integrand(min, integrandBuffer);
			resultBuffer.multAdd(*umin - min, integrandBuffer);

			integrand(*umin, integrandBuffer);
			resultBuffer.multAdd(*(umin + 1) - min, integrandBuffer);

			while (++umin != max)
			{
				integrand(*umin, integrandBuffer);
				resultBuffer.multAdd(*(umin + 1) - *(umin - 1), integrandBuffer);
			}

			integrand(*umin, integrandBuffer);
			resultBuffer.multAdd(*umin - *(umin - 1), integrandBuffer);

			resultBuffer *= 0.5f;
		}
		else
		{
			integrand(min, integrandBuffer);
			resultBuffer += integrandBuffer;

			integrand(*umin, integrandBuffer);
			resultBuffer += integrandBuffer;

			resultBuffer *= 0.5f * (*umin - min);
		}
	}

	/**
	 * @brief One-dimensional trapezoidal integration in one-dimensional frequency space. 
	 * First argument of the integrand function is the frequency value, the second argument is the return value of the integrand, passed by reference. 
	 * 
	 * @tparam T Integrand type. 
	 * @param[in] min Iterator to lower boundary frequency value. 
	 * @param[in] max Upper boundary frequency value. 
	 * @param[in] integrand Integrand function. 
	 * @param[out] integrandBuffer Return value buffer of the integrand. 
	 * @param[out] resultBuffer Value of the integral. 
	 */
	template <class T>
	void integrateWithObscureRightBoundary(const FrequencyIterator min, const float max, const std::function<void(float, T &)> &integrand, T &integrandBuffer, T &resultBuffer)
	{
		ASSERT(&integrandBuffer != &resultBuffer);
		ASSERT(*min <= max, "Lower integration boundary must not be larger than upper boundary. ");

		resultBuffer.reset();
		FrequencyIterator umin(min);
		FrequencyIterator umax = FrgCommon::frequency().lesser(max);
		if (umin != umax)
		{
			integrand(*umin, integrandBuffer);
			resultBuffer.multAdd(*(umin + 1) - *umin, integrandBuffer);

			while (++umin != umax)
			{
				integrand(*umin, integrandBuffer);
				resultBuffer.multAdd(*(umin + 1) - *(umin - 1), integrandBuffer);
			}

			integrand(*umin, integrandBuffer);
			resultBuffer.multAdd(max - *(umin - 1), integrandBuffer);

			integrand(max, integrandBuffer);
			resultBuffer.multAdd(max - *umin, integrandBuffer);

			resultBuffer *= 0.5f;
		}
		else
		{
			integrand(max, integrandBuffer);
			resultBuffer += integrandBuffer;

			integrand(*umin, integrandBuffer);
			resultBuffer += integrandBuffer;

			resultBuffer *= 0.5f * (max - *umin);
		}
	}

	/**
	 * @brief One-dimensional trapezoidal integration in one-dimensional frequency space. 
	 * First argument of the integrand function is the frequency value, the second argument is the return value of the integrand, passed by reference. 
	 * 
	 * @tparam T Integrand type. 
	 * @param[in] min Lower boundary frequency value. 
	 * @param[in] max Upper boundary frequency value. 
	 * @param[in] integrand Integrand function. 
	 * @param[out] integrandBuffer Return value buffer of the integrand. 
	 * @param[out] resultBuffer Value of the integral. 
	 */
	template <class T>
	void integrateWithObscureBoundaries(const float min, const float max, const std::function<void(float, T &)> &integrand, T &integrandBuffer, T &resultBuffer)
	{
		ASSERT(&integrandBuffer != &resultBuffer);
		ASSERT(min <= max, "Lower integration boundary must not be larger than upper boundary. ");

		resultBuffer.reset();
		FrequencyIterator umin = FrgCommon::frequency().greater(min);
		FrequencyIterator umax = FrgCommon::frequency().lesser(max);
		if (umax >= umin)
		{
			integrand(min, integrandBuffer);
			resultBuffer.multAdd(*umin - min, integrandBuffer);

			if (umax != umin)
			{
				integrand(*umin, integrandBuffer);
				resultBuffer.multAdd(*(umin + 1) - min, integrandBuffer);

				while (++umin != umax)
				{
					integrand(*umin, integrandBuffer);
					resultBuffer.multAdd(*(umin + 1) - *(umin - 1), integrandBuffer);
				}

				integrand(*umin, integrandBuffer);
				resultBuffer.multAdd(max - *(umin - 1), integrandBuffer);
			}
			else
			{
				integrand(*umin, integrandBuffer);
				resultBuffer.multAdd(max - min, integrandBuffer);
			}

			integrand(max, integrandBuffer);
			resultBuffer.multAdd(max - *umin, integrandBuffer);

			resultBuffer *= 0.5f;
		}
		else
		{
			integrand(max, integrandBuffer);
			resultBuffer += integrandBuffer;

			integrand(min, integrandBuffer);
			resultBuffer += integrandBuffer;

			resultBuffer *= 0.5f * (max - min);
		}
	}
}