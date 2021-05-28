/**
 * @file LoadManager.hpp
 * @author Finn Lasse Buessen
 * @brief Implementation of an automatic hybrid OpenMP / MPI load balancer. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <vector>
#include <functional>
#include <thread>
#include <mutex>
#include <boost/date_time.hpp>
#include "lib/Log.hpp"
#include "lib/Exception.hpp"

#ifndef DISABLE_MPI
#include "mpi.h"
#endif
#ifndef DISABLE_OMP
#include "omp.h"
#endif

#define HMP_CHUNK_PROPERTY_STACK 0 ///< Memory offset of the stack id in the chunk properties. 
#define HMP_CHUNK_PROPERTY_BEGIN 1 ///< Memory offset of the workload begin in the chunk properties. 
#define HMP_CHUNK_PROPERTY_END 2 ///< Memory offset of the workload end in the chunk properties. 

#ifndef DISABLE_MPI
#define HMP_MPI_ENABLED ///< Defined, if MPI parallelization is enabled. 
#define HMP_ENABLE_IF_MPI(X) X ///< Print argument if MPI parallelization is enabled. 
#define HMP_DISABLE_IF_MPI(X) ///< Do not print argument if MPI parallelization is enabled. 
#define HMP_APPEND_IF_MPI(X) ,X ///< Append argument if MPI parallelizatino is enabled. 
#else
#define HMP_ENABLE_IF_MPI(X) ///< Print argument if MPI parallelization is enabled. 
#define HMP_DISABLE_IF_MPI(X) X ///< Do not print argument if MPI parallelization is enabled. 
#define HMP_APPEND_IF_MPI(X) ///< Append argument if MPI parallelizatino is enabled. 
#endif

namespace HMP
{
	class LoadManager;
	typedef int StackIdentifier; ///< DataStack identifier. 

	/**
	 * @brief Common load manager interface for both, the MPI server rank and slave ranks. 
	 * @details This class provides an implementation of an automatic load balancing scheme for hybrid shared/distributed memory parallelization. 
	 * 	Intra-node parallelization is implemented via OpenMP. Inter-node parallelization is implemented via MPI. 
	 *  A LoadManager instance is created on every MPI rank by calling the LoadManager::newLoadManager routine; 
	 *  This will create a LoadManagerMaster instance on the designated MPI master rank, and LoadManagerSlave instances on all remaining ranks. 
	 * 	The work load is registered with the LoadManager in the form of DataStacks. Each data stack is associated with an array of numbers 
	 *  and with a functional description (referred to as calculators) on how to perform computations on the individual array entries. 
	 *  Different types of data stack exist with different calculator properties, see DataStackBase::StackType for more information. 
	 *  The LoadManager interface provides functions to create and register the different stack types with the LoadManager. 
	 *  
	 *  The role of the LoadManagerMaster is to break down the collective computation of all array entries into smaller chunks of work, 
	 *  to distribute them either locally or across multiple MPI ranks, and to retrieve the result on the master rank. 
	 *  In order to execute the calculators on a specific stack, the LoadManager provides the LoadManager::calculate interface. 
	 *  In order to make the results available also on other MPI ranks, the LoadManager::broadcast interface is provided. 
	 *  Various different types of data stacks exist, see the appropriate member functions for adding new stacks. 
	 * 
	 * @see LoadManager::newLoadManager
	 * @see LoadManager::calculate
	 * @see LoadManager::broadcast
	 */
	class LoadManager
	{
	protected:
		/**
		 * @brief Abstract base class for DataStack types. 
		 * @details This abstract type does not actually store data or prescriptions for calculations. 
		 * It merely specifies the virtual interface for inheritance. Subclasses should then implement data storage and calculators, see concrete implementations for more details. 
		 * 
		 * @see DataStack
		 */
		struct DataStackBase
		{
			/**
			 * @brief List to identify the different traits of DataStack. 
			 */
			enum struct StackType
			{
				/**
				 * @brief Explicit DataStack 
				 * @details An explicit DataStack contains an explicit calculator, i.e. a function with a single integer as an argument, and a non-void return type. 
				 * When the calculator is invoked, its return value is written to the data element whose position is specified by the integer argument. 
				 * @see LoadManager::addMasterStackExplicit
				 */
				Explicit, 
				/**
				 * @brief Implicit DataStack 
				 * @details An implicit DataStack contains an implicit calculator, i.e. a function with a single integer as an argument, and a void return type. 
				 * When the calculator is invoked, it is passed the index of a data element stored in the stack. The calculator is expected to modify that element. 
				 * @see LoadManager::addMasterStackImplicit
				 */
				Implicit,
				/**
				 * @brief Slave DataStack 
				 * @details A slave DataStack does not contain a calculator. Instead, it is associated with an implicit data stack. 
				 * The data storage is expected to have the same size as the associated implicit stack, and the associated stack's implicit calculator is expected to modify also the slave stack's element with the corresponding index. 
				 * Whenever data of the associated implicit stack is communicated between different MPI ranks, the corresponding slave stack data is also communicated. 
				 * @see LoadManager::addSlaveStack
				 */
				Slave,
				/**
				 * @brief Passive DataStack 
				 * @details A passive DataStack does not contain a calculator. Changes to its data can only be made on the MPI server rank and need to be manually broadcasted. 
				 * @see LoadManager::addPassiveStack
				 */
				Passive
			};

			/**
			 * @brief Internal message tags for MPI communication. 
			 */
			enum struct MessageTag
			{
				Chunk, ///< MPI message contains Chunk information
				ChunkReturn, ///< MPI message contains the result of a chunk computation
			};

			/**
			 * @brief Virtual destructor.
			 */
			virtual ~DataStackBase() {};

			/**
			 * @brief Virtual function to apply an internal calculator to the n-th value stored in the stack. 
			 * 
			 * @param n Specifies the index to which the calculator should be applied. 
			 */
			virtual void applyCalculator(const int n) {};

			#ifdef HMP_MPI_ENABLED
			/**
			 * @brief Virtual function to send a data block to a DataStack on a different MPI rank, typically the server rank. 
			 * 
			 * @param offset Offset to the beginning of the data block to be sent, measured in number of entries (or number of entry tuples, if DataStackBase::typeMultiplicity is greater than one).
			 * @param count Number of entries to be sent. 
			 * @param serverRank Receiver's MPI rank, typically the server rank. 
			 * @param communicator The MPI communicator used for communication. 
			 */
			virtual void send(const int offset, const int count, const int serverRank, const MPI_Comm communicator) const {};

			/**
			 * @brief Virtual function to asynchronously receive a data block from a different MPI rank. 
			 * 
			 * @param[in] offset Offset to the beginning of the data block to be received, measured in number of entries (or number of entry tuples, if DataStackBase::typeMultiplicity is greater than one).
			 * @param[in] count Number of entries to be received. 
			 * @param[in] rank Sender's MPI rank. 
			 * @param[in] communicator The MPI communicator used for communication. 
			 * @param[out] request MPI request object for the communication; Should be used to determine whether the non-blocking receive has been completed. 
			 */
			virtual void receive(const int offset, const int count, const int rank, const MPI_Comm communicator, MPI_Request &request) {};

			/**
			 * @brief Virtual function to broadcast a data block from the server rank to all other MPI ranks. 
			 * 
			 * @param serverRank The MPI rank to take the role of the sender. 
			 * @param communicator The MPI communicator used for communication. 
			 */
			virtual void broadcast(const int serverRank, const MPI_Comm communicator) {};
			#endif

			StackType type; ///< Specifies the trait of the stack @see StackType.
			StackIdentifier master; ///< Specifies the id of an associated stack. Only relevant for slave stacks, otherwise set to -1. @see StackType::Slave
			int size; ///< Number of elements (or tuples of elements if DataStackBase::typeMultiplicity is greater than one) stored in the stack. 
			int typeMultiplicity; ///< Multiplicity of each element. Cannot be greater than one for explicit stacks. If greater than one, the data stack is assumed to consist of tuples of fundamental data types. Element indexing then refers to the tuples, not the fundamental data types. 
			int recommendedChunkSizeMultiple; ///< When breaking the data stack down into smaller work chunks, attempt to form chunks whose size is a multiple of the given value. This is helpful if calculators vary in runtime, but can be joined to groups whose collective runtime is expected to be constant. 
			int recommendedChunksPerRank; ///< When breaking the data stack down into smaller work chunks, attempt to form approximately the specified number of chunks per MPI rank. 
			bool autoBroadcast; ///< If set to true, modifications to the stack's data that are a consequence of the onvication of calculators are automatically communicated across all MPI ranks. If set to false, they are only sent to the MPI server rank. 
		};

		/**
		 * @brief Concrete derivation of DataStackBase, managing an array of fundamental data types (or tuples thereof, if DataStackBase::typeMultiplicity is greater than one). 
		 * 
		 * @tparam StackT The fundamental data type stored in the internal data array. 
		 */
		template <class StackT> struct DataStack : public DataStackBase
		{
			friend class LoadManager;
		protected:
			/**
			 * @brief Construct a new Data Stack object. 
			 * 
			 * @see LoadManager::addMasterStackExplicit
			 * @see LoadManager::addMasterStackImplicit
			 * @see LoadManager::addSlaveStack
			 * @see LoadManager::addPassiveStack
			 */
			DataStack() {}

			/**
			 * @brief Invoke internal calculator and write result to the data array at the specified index. 
			 * 
			 * @param i Index of the element to which the result is written. 
			 */
			void applyCalculator(const int i) override
			{
				if (type == StackType::Implicit) implicitCalculator(i);
				else if (type == StackType::Explicit) data[i] = explicitCalculator(i);
			}

			#ifdef HMP_MPI_ENABLED
			/**
			 * @brief Send a data block to a different MPI rank, typically the server rank. 
			 * 
			 * @param offset Id of the first element of the data block to be sent, measured in number of entries (or number of entry tuples, if DataStackBase::typeMultiplicity is greater than one).
			 * @param count Number of elements (or element tuples) to be sent. 
			 * @param serverRank Receiver's MPI rank, typically the server rank. 
			 * @param communicator The MPI communicator used for communication. 
			 */
			void send(const int offset, const int count, const int serverRank, const MPI_Comm communicator) const override
			{
				MPI_Send(static_cast<void *>(data + typeMultiplicity * offset), typeMultiplicity *count * sizeof(StackT), MPI_BYTE, serverRank, static_cast<int>(MessageTag::ChunkReturn), communicator);
			}

			/**
			 * @brief Asynchronously receive a data block from a different MPI rank. 
			 * 
			 * @param[in] offset Id of the first element of the data block to be received, measured in number of entries (or number of entry tuples, if DataStackBase::typeMultiplicity is greater than one).
			 * @param[in] count Number of elements (or element tuples) to be received. 
			 * @param[in] rank Sender's MPI rank. 
			 * @param[in] communicator The MPI communicator used for communication. 
			 * @param[out] request MPI request object for the communication; Should be used to determine whether the non-blocking receive has been completed. 
			 */
			void receive(const int offset, const int count, const int rank, const MPI_Comm communicator, MPI_Request &request) override
			{
				MPI_Irecv(static_cast<void *>(data + typeMultiplicity * offset), typeMultiplicity *count * sizeof(StackT), MPI_BYTE, rank, static_cast<int>(MessageTag::ChunkReturn), communicator, &request);
			}

			/**
			 * @brief Broadcast a data block from the server rank to all other MPI ranks. 
			 * 
			 * @param serverRank The MPI rank to take the role of the sender. 
			 * @param communicator The MPI communicator used for communication. 
			 */
			void broadcast(const int serverRank, const MPI_Comm communicator) override
			{
				MPI_Bcast(static_cast<void *>(data), typeMultiplicity * size * sizeof(StackT), MPI_BYTE, serverRank, communicator);
			}
			#endif

			std::function<StackT(int)> explicitCalculator; ///< Explicit calculator. Only relevant if DataStackBase::type is set to StackType::Explicit. 
			std::function<void(int)> implicitCalculator; ///< Implicit calculator. Only relevant if DataStackBase::type is set to StackType::Implicit. 
			StackT *data; ///< Internal data array. 
		};

		/**
		 * @brief Definition of a small chunk of work. Used to communicate workload between different LoadManager instances. 
		 */
		struct Chunk
		{
		public:
			///Construct an empty work chunk
			Chunk()
			{
				properties[HMP_CHUNK_PROPERTY_STACK] = -1;
				properties[HMP_CHUNK_PROPERTY_BEGIN] = 0;
				properties[HMP_CHUNK_PROPERTY_END] = 0;
			};

			/**
			 * @brief Construct a work chunk for a specific workload. 
			 * 
			 * @param stackId Id of the stack that the work should be performed on. 
			 * @param begin Index of the first element whose calculator should be invoked. 
			 * @param end Index of the last element whose calculator should be invoked. 
			 */
			Chunk(const int stackId, const int begin, const int end)
			{
				properties[HMP_CHUNK_PROPERTY_STACK] = stackId;
				properties[HMP_CHUNK_PROPERTY_BEGIN] = begin;
				properties[HMP_CHUNK_PROPERTY_END] = end;
			}

			///Check of the workload is empty
			bool isVoid() const
			{
				return (properties[HMP_CHUNK_PROPERTY_STACK] == -1);
			}

			///Comparison operator
			bool operator==(const Chunk &rhs) const
			{
				return (properties[HMP_CHUNK_PROPERTY_STACK] == rhs.properties[HMP_CHUNK_PROPERTY_STACK]) && (properties[HMP_CHUNK_PROPERTY_BEGIN] == rhs.properties[HMP_CHUNK_PROPERTY_BEGIN]) && (properties[HMP_CHUNK_PROPERTY_END] == rhs.properties[HMP_CHUNK_PROPERTY_END]);
			}

			///Negative comparison operator
			bool operator!=(const Chunk &rhs) const
			{
				return !this->operator==(rhs);
			}

			int properties[3]; ///< Workload specification. First value describes the stack id, second value describes the first element of the workload, and the third value the last element of the workload.  
		};

	public:
		///Destroy the LoadManager object
		virtual ~LoadManager()
		{
			HMP_ENABLE_IF_MPI(MPI_Comm_free(&_communicator));

			while (_stacks.size() > 0)
			{
				delete _stacks.back();
				_stacks.pop_back();
			}
		}

		/**
		 * @brief Create an explicit data stack and attach it to the LoadManager. 
		 * 
		 * @tparam StackT The fundamental data type of the DataStack to be created. 
		 * @param data Data array on which the stack operates. The allocated size should be at least size * typeof(StackT). 
		 * @param size Number of elements in the stack.
		 * @param calculator Calculator object associated with the stack. 
		 * @param recommendedChunkSizeMultiple Suggested quantization of elements per work chunk. 
		 * @param recommendedChunksPerRank Suggested number of work chunks per MPI rank. 
		 * @param autoBroadcast Enable or disable auto broadcast. 
		 * @return StackIdentifier Id of the newly generated stack as registered with the LoadManager. 
		 * 
		 * @see DataStackBase::StackType::Explicit
		 * @see DataStack
		 */
		template <class StackT> StackIdentifier addMasterStackExplicit(StackT *const data, const int size, const std::function<StackT(int)> &calculator, const int recommendedChunkSizeMultiple = 1, const int recommendedChunksPerRank = 10, const bool autoBroadcast = false)
		{
			DataStack<StackT> *ds = new DataStack<StackT>();
			ds->type = DataStackBase::StackType::Explicit;
			ds->master = -1;
			ds->size = size;
			ds->typeMultiplicity = 1;
			ds->recommendedChunkSizeMultiple = recommendedChunkSizeMultiple;
			ds->recommendedChunksPerRank = recommendedChunksPerRank;
			ds->autoBroadcast = autoBroadcast;
			ds->explicitCalculator = calculator;
			ds->data = data;
			return _registerStack(ds);
		}

		/**
		 * @brief Create an implicit data stack and attach it to the LoadManager. 
		 * 
		 * @tparam StackT The fundamental data type of the DataStack to be created. 
		 * @param data Data array on which the stack operates. The allocated size should be at least size * typeMultiplicity * typeof(StackT). 
		 * @param size Number of elements (or element tuples) in the stack.
		 * @param calculator Calculator object associated with the stack. 
		 * @param typeMultiplicity Element tuple size. 
		 * @param recommendedChunkSizeMultiple Suggested quantization of elements per work chunk. 
		 * @param recommendedChunksPerRank Suggested number of work chunks per MPI rank. 
		 * @param autoBroadcast Enable or disable auto broadcast. 
		 * @return StackIdentifier Id of the newly generated stack as registered with the LoadManager. 
		 * 
		 * @see DataStackBase::StackType::Implicit
		 * @see DataStack
		 */
		template <class StackT> StackIdentifier addMasterStackImplicit(StackT *const data, const int size, const std::function<void(int)> &calculator, const int typeMultiplicity = 1, const int recommendedChunkSizeMultiple = 1, const int recommendedChunksPerRank = 10, const bool autoBroadcast = false)
		{
			DataStack<StackT> *ds = new DataStack<StackT>();
			ds->type = DataStackBase::StackType::Implicit;
			ds->master = -1;
			ds->size = size;
			ds->typeMultiplicity = typeMultiplicity;
			ds->recommendedChunkSizeMultiple = recommendedChunkSizeMultiple;
			ds->recommendedChunksPerRank = recommendedChunksPerRank;
			ds->autoBroadcast = autoBroadcast;
			ds->implicitCalculator = calculator;
			ds->data = data;
			return _registerStack(ds);
		}

		/**
		 * @brief Create a slave data stack and attach it to the LoadManager. 
		 * 
		 * @tparam StackT The fundamental data type of the DataStack to be created. 
		 * @param data Data array on which the stack operates. The allocated size should be at least size * typeMultiplicity * typeof(StackT). 
		 * @param size Number of elements (or element tuples) in the stack.
		 * @param master Stack id of the associated implicit master stack. 
		 * @param typeMultiplicity Element tuple size. 
		 * @return StackIdentifier Id of the newly generated stack as registered with the LoadManager. 
		 * 
		 * @see DataStackBase::StackType::Slave
		 * @see DataStack
		 */
		template <class StackT> StackIdentifier addSlaveStack(StackT *const data, const int size, const StackIdentifier master, const int typeMultiplicity = 1)
		{
			DataStack<StackT> *ds = new DataStack<StackT>();
			ds->type = DataStackBase::StackType::Slave;
			ds->master = master;
			ds->size = size;
			ds->typeMultiplicity = typeMultiplicity;
			ds->autoBroadcast = false;
			ds->data = data;
			return _registerStack(ds);
		}

		/**
		 * @brief Create a passive data stack and attach it to the LoadManager. 
		 * 
		 * @tparam StackT The fundamental data type of the DataStack to be created. 
		 * @param data Data array on which the stack operates. The allocated size should be at least size * typeMultiplicity * typeof(StackT). 
		 * @param size Number of elements (or element tuples) in the stack.
		 * @return StackIdentifier Id of the newly generated stack as registered with the LoadManager. 
		 * 
		 * @see DataStackBase::StackType::Passive
		 * @see DataStack
		 */
		template <class StackT> StackIdentifier addPassiveStack(StackT *const data, const int size)
		{
			DataStack<StackT> *ds = new DataStack<StackT>();
			ds->type = DataStackBase::StackType::Passive;
			ds->master = -1;
			ds->size = size;
			ds->typeMultiplicity = 1;
			ds->autoBroadcast = false;
			ds->data = data;
			return _registerStack(ds);
		}

		/**
		 * @brief Calculate a list of stacks, where the stack identifiers are provided in list form. 
		 * 
		 * @param stackIds Pointer to the first StackIdentifier. 
		 * @param size Number of stacks. 
		 */
		virtual void calculate(const StackIdentifier *stackIds, const int size) = 0;

		/**
		 * @brief Calculate a list of stacks, where the stack identifiers are provided in initializer list form. 
		 * 
		 * @param stackIds Initializer list of StackIdentifiers.
		 */
		void calculate(const std::initializer_list<StackIdentifier> &stackIds)
		{
			calculate(stackIds.begin(), int(stackIds.size()));
		}

		/**
		 * @brief Calculate a single stack.
		 * 
		 * @param stackId StackIdentifier of the stack to be calculated. 
		 */
		void calculate(const StackIdentifier stackId)
		{
			calculate({ stackId });
		}

		/**
		 * @brief Calculate all stacks. 
		 */
		void calculateAll()
		{
			std::vector<StackIdentifier> all(_stacks.size());
			for (StackIdentifier i = 0; i < StackIdentifier(_stacks.size()); ++i) all[i] = i;
			calculate(all.data(), int(all.size()));
		}

		/**
		 * @brief Broadcast a list of stacks, where the stack identifiers are provided in list form. 
		 * 
		 * @param stackIds Pointer to the first StackIdentifier. 
		 * @param size Number of stacks. 
		 */
		void broadcast(const StackIdentifier *stackIds, const int size)
		{
			#ifdef HMP_MPI_ENABLED
			for (int i = 0; i < size; ++i)
			{
				for (StackIdentifier s = 0; s < StackIdentifier(_stacks.size()); ++s)
				{
					if (s == stackIds[i] || _stacks[s]->master == stackIds[i]) _stacks[s]->broadcast(_serverRank, _communicator);
				}
			}
			#endif
		}

		/**
		 * @brief Broadcast a list of stacks, where the stack identifiers are provided in initializer list form. 
		 * 
		 * @param stackIds Initializer list of StackIdentifiers.
		 */
		void broadcast(const std::initializer_list<StackIdentifier> &stackIds)
		{
			broadcast(stackIds.begin(), int(stackIds.size()));
		}

		/**
		 * @brief Broadcast a single stack. 
		 * 
		 * @param stackId StackIdentifier of the stack to be broadcasted. 
		 */
		void broadcast(const StackIdentifier stackId)
		{
			broadcast({ stackId });
		}

		/**
		 * @brief Broadcast all stacks. 
		 */
		void broadcastAll()
		{
			std::vector<StackIdentifier> all(_stacks.size());
			for (StackIdentifier i = 0; i < StackIdentifier(_stacks.size()); ++i) all[i] = i;
			broadcast(all.data(), int(all.size()));
		}

		/**
		 * @brief Virtual implementation to print runtime statistics, including information on the efficiency of LoadManager instances running on different MPI ranks. 
		 */
		virtual void printRuntimeStatistics() const {}

	protected:
		/**
		 * @brief Construct a new Load Manager object
		 * 
		 * @param serverRank The MPI rank to take on the master role. 
		 * @param communicator The MPI communicator to operate on. 
		 */
		LoadManager(const int serverRank HMP_APPEND_IF_MPI(const MPI_Comm communicator)) 
		{
			HMP_ENABLE_IF_MPI(MPI_Comm_dup(communicator, &_communicator));
			HMP_ENABLE_IF_MPI(MPI_Comm_size(_communicator, &_commSize));
			HMP_DISABLE_IF_MPI(_commSize = 1);
			HMP_ENABLE_IF_MPI(MPI_Comm_rank(_communicator, &_rank));
			HMP_DISABLE_IF_MPI(_rank = 0);

			if (serverRank >= _commSize) throw Exception(Exception::Type::MpiError, "Server rank must not exceed communicator size.");
			else _serverRank = serverRank;
		}

		/**
		 * @brief Register a DataStackBase with the LoadManager. By registering the stack, the LoadManager assumes responsibility to synchronize data between MPI ranks as necessary. 
		 * 
		 * @param stack The stack to be registered. 
		 * @return StackIdentifier Unique identifier of the newly registered stack; Used e.g. for LoadManager::calculate calls. 
		 */
		virtual StackIdentifier _registerStack(DataStackBase *stack)
		{
			//If stack is a slave stack, check that the type multiplicity matches the type multiplicity of the master
			if (stack->type == DataStackBase::StackType::Slave && stack->typeMultiplicity != _stacks[stack->master]->typeMultiplicity) throw Exception(Exception::Type::ArgumentError, "Stack type multiplicity of slave stack does not match type multiplicity of associated master stack.");

			_stacks.push_back(stack);
			return StackIdentifier(_stacks.size() - 1);
		}

		/**
		 * @brief Calculate the workload defined by a specific chunk. 
		 * 
		 * @param chunk Provides the workload definition. 
		 */
		void _calculateChunk(const Chunk &chunk)
		{
			#ifndef DISABLE_OMP
			#pragma omp parallel for schedule(guided)
			#endif
			for (int i = chunk.properties[HMP_CHUNK_PROPERTY_BEGIN]; i < chunk.properties[HMP_CHUNK_PROPERTY_END]; ++i) _stacks[chunk.properties[HMP_CHUNK_PROPERTY_STACK]]->applyCalculator(i);
		}

		std::vector<DataStackBase *> _stacks; ///< List of all registered stacks. 
		int _serverRank; ///< Designated MPI master rank. 
		int _rank; ///< MPI rank of the current LoadManager instance. 
		int _commSize; ///< MPI communicator size. 
		HMP_ENABLE_IF_MPI(MPI_Comm _communicator); ///< MPI communicator to operate on. 
	};

	/**
	 * @brief LoadManager master implementation, which is responsible for distributing work and synchronizing data. 
	 * 
	 * @see LoadManager
	 */
	class LoadManagerMaster : public LoadManager
	{
		friend LoadManager *newLoadManager(const int serverRank HMP_APPEND_IF_MPI(const MPI_Comm communicator));
	public:
		/**
		 * @brief Calculate a list of stacks, where the stack identifiers are provided in list form. 
		 * 
		 * @param stackIds Pointer to the first StackIdentifier. 
		 * @param size Number of stacks. 
		 */
		virtual void calculate(const StackIdentifier *stackIds, const int size) override
		{
			//reset calculation runtime statistics
			for (int i = 0; i < _commSize; ++i)
			{
				for (StackIdentifier s = 0; s < StackIdentifier(_stacks.size()); ++s)
				{
					_currentCalculationComputeTimeBuffer[i * _stacks.size() + s] = 0.0f;
					_currentCalculationWorkDone[i][s] = 0;
					_currentCalculationTime[i][s] = 0.0f;

				}
			}
			boost::posix_time::ptime tic = boost::posix_time::microsec_clock::local_time();

			//init chunk spawner
			_initChunkSpawner(stackIds, size);

			//spawn local worker
			std::thread *t = new std::thread([&]() { this->_runLocalClient(); });

			#ifdef HMP_MPI_ENABLED
			//issue initial chunks
			for (int i = 0; i < _commSize; ++i) if (i != _serverRank) _issueChunk(i);

			//collect results and issue consecutive chunks
			int receiveFlag;
			MPI_Status receiveStatus;
			for (;;)
			{
				begin:
				for (int rank = 0; rank < _commSize; ++rank)
				{
					if (rank == _serverRank || _pendingRequests[rank] == MPI_REQUEST_NULL) continue;
					MPI_Test(_pendingRequests + rank, &receiveFlag, &receiveStatus);
					if (receiveFlag == 1)
					{
						_despawnChunk(rank);
						_issueChunk(receiveStatus.MPI_SOURCE);
					}
				}
				for (int rank = 0; rank < _commSize; ++rank) if (_pendingRequests[rank] != MPI_REQUEST_NULL) goto begin;
				break;
			}
			#endif

			//join local worker
			t->join();

			//broadcast result
			for (int s = 0; s < size; ++s)
			{
				if (_stacks[stackIds[s]]->autoBroadcast) broadcast(stackIds[s]);
			}

			//gather runtime statistics
			HMP_ENABLE_IF_MPI(MPI_Gather(MPI_IN_PLACE, int(_stacks.size()), MPI_FLOAT, _currentCalculationComputeTimeBuffer.data(), int(_stacks.size()), MPI_FLOAT, _serverRank, _communicator));

			for (int i = 0; i < _commSize; ++i)
			{
				for (StackIdentifier s = 0; s < StackIdentifier(_stacks.size()); ++s)
				{
					_totalComputeTime[i][s] += _currentCalculationComputeTimeBuffer[i * _stacks.size() + s];
				}
			}

			boost::posix_time::ptime toc = boost::posix_time::microsec_clock::local_time();
			_totalCalculationTime += float((toc - tic).total_milliseconds());
		}

		/**
		 * @brief Print runtime statistics, including information on the efficiency of LoadManager instances running on different MPI ranks. 
		 */
		void printRuntimeStatistics() const override
		{
			Log::log << Log::LogLevel::Debug << "LoadManager total time spent in calculate() calls is " << _totalCalculationTime << "ms" << Log::endl;
			for (int i = 0; i < _commSize; ++i)
			{
				float timePerRank = 0.0;
				for (StackIdentifier j = 0; j < StackIdentifier(_stacks.size()); ++j) timePerRank += _totalComputeTime[i][j];
				Log::log << Log::LogLevel::Debug << "LoadManager (rank " << i << ") active computing time was " << std::setiosflags(std::ios::fixed) << std::setprecision(0) << timePerRank << "ms (Efficiency: " << std::setprecision(1) << 100.0f * timePerRank / _totalCalculationTime << "%)" << Log::endl;

				for (StackIdentifier j = 0; j < StackIdentifier(_stacks.size()); ++j) Log::log << Log::LogLevel::Debug << "\t active computing time on stack " << j << " was " << _totalComputeTime[i][j] << "ms" << Log::endl;
			}
		}

	protected:
		/**
		 * @brief Construct a new LoadManagerMaster object
		 * 
		 * @param serverRank The MPI rank to take on the master role. 
		 * @param communicator The MPI communicator used for communication. 
		 */
		LoadManagerMaster(const int serverRank HMP_APPEND_IF_MPI(const MPI_Comm communicator)) : LoadManager(serverRank HMP_APPEND_IF_MPI(communicator))
		{
			_totalCalculationTime = 0.0f;
			_totalComputeTime = new std::vector<float>[_commSize];
			_currentCalculationWorkDone = new std::vector<int>[_commSize];
			_currentCalculationTime = new std::vector<float>[_commSize];
			_currentCalculationChunkSpawntime = new boost::posix_time::ptime[_commSize];
			_currentCalculationChunkSpawned = new Chunk[_commSize];

			#ifdef HMP_MPI_ENABLED
			_pendingRequests = new MPI_Request[_commSize];
			for (int i = 0; i < _commSize; ++i) _pendingRequests[i] = MPI_REQUEST_NULL;
			#endif
		}

		/**
		 * @brief Destroy the LoadManagerMaster object
		 */
		~LoadManagerMaster()
		{
			delete[] _totalComputeTime;
			delete[] _currentCalculationWorkDone;
			delete[] _currentCalculationTime;
			delete[] _currentCalculationChunkSpawntime;
			delete[] _currentCalculationChunkSpawned;
			
			HMP_ENABLE_IF_MPI(delete[] _pendingRequests);
		}

		/**
		 * @brief Register a DataStackBase with the LoadManager. By registering the stack, the LoadManager assumes responsibility to synchronize data between MPI ranks as necessary. 
		 * 
		 * @param stack The stack to be registered. 
		 * @return StackIdentifier Unique identifier of the newly registered stack; Used e.g. for LoadManager::calculate calls. 
		 */
		virtual StackIdentifier _registerStack(DataStackBase *stack) override
		{
			StackIdentifier identifier = LoadManager::_registerStack(stack);

			for (int i = 0; i < _commSize; ++i)
			{
				_totalComputeTime[i].push_back(0.0f);
				_currentCalculationWorkDone[i].push_back(0);
				_currentCalculationTime[i].push_back(0.0f);
			}

			_currentCalculationComputeTimeBuffer.resize(_commSize * _stacks.size());
			_currentCalculationStackMask.resize(_stacks.size());
			_currentCalculationStackProgress.resize(_stacks.size());

			return identifier;
		}

		/**
		 * @brief Worker loop to run calculators locally. 
		 */
		void _runLocalClient()
		{
			for (;;)
			{
				//get chunk
				Chunk c = _spawnChunk(_serverRank);
				if (c.isVoid()) break;

				boost::posix_time::ptime tic = boost::posix_time::microsec_clock::local_time();
				_calculateChunk(c);
				boost::posix_time::ptime toc = boost::posix_time::microsec_clock::local_time();
				_currentCalculationComputeTimeBuffer[c.properties[HMP_CHUNK_PROPERTY_STACK]] += (toc - tic).total_milliseconds();

				_despawnChunk(_serverRank);
			}
		}

		/**
		 * @brief Prepare the chunk spawner, i.e. the routine responsible for breaking down workload into smaller chunsk, for the execution of calculators on the specified list of stacks. 
		 * 
		 * @param stackIds Pointer to the first StackIdentifier to be calculated. 
		 * @param size Number of StackIdentifiers to be calculated. 
		 */
		void _initChunkSpawner(const StackIdentifier *stackIds, const int size)
		{
			//reset calculation statistics
			for (int i = 0; i < _commSize; ++i)
			{
				for (StackIdentifier s = 0; s < StackIdentifier(_stacks.size()); ++s)
				{
					_currentCalculationWorkDone[i][s] = 0;
					_currentCalculationTime[i][s] = 0.0f;
				}
			}
			for (StackIdentifier s = 0; s < StackIdentifier(_stacks.size()); ++s)
			{
				_currentCalculationStackMask[s] = false;
				_currentCalculationStackProgress[s] = 0;
			}

			//enable selected stacks
			for (int i = 0; i < size; ++i) _currentCalculationStackMask[stackIds[i]] = true;
		}

		/**
		 * @brief Prepare the chunk spawner, i.e. the routine responsible for breaking down workload into smaller chunsk, for the execution of calculators on the specified initializer list of stacks.
		 * 
		 * @param stackIds Initializer list of the stacks to be calculated. 
		 */
		void _initChunkSpawner(const std::initializer_list<StackIdentifier> &stackIds)
		{
			_initChunkSpawner(stackIds.begin(), int(stackIds.size()));
		}

		/**
		 * @brief The chunk spawner routine, which is responsible for breaking down the overall workload of calculator executions into smaller chunks of workload. 
		 *  The algorithm to spawn a chunk and determine its workload is as follows:
		 *  1. Choose a stack that has not finished computing. 
		 *  2. Determine the maximum chunk size according to the stack's recommended number of chunks per rank. 
		 *  3. Round up that number to be a multiple of the recommended chunk size. 
		 *  4. Determine chunk size according to relative computing power of the different ranks, based on performance on previous chunks. 
		 *  5. Clip chunk size to minimum of 100ms expected return time and maximum as determined before. 
		 * 
		 * @param rank MPI rank for which the workload chunk is being requested. 
		 * @return Chunk Definition of the workload. 
		 */
		Chunk _spawnChunk(const int rank)
		{
			std::lock_guard<std::mutex> lock(_currentCalculationChunkSpawnerLock);

			Chunk c;
			for (StackIdentifier s = 0; s < StackIdentifier(_stacks.size()); ++s)
			{
				//skip inactive and completed stacks
				if (!_currentCalculationStackMask[s] || _currentCalculationStackProgress[s] >= _stacks[s]->size) continue;

				//determine max chunk size
				int maximumWorkShare = int(_stacks[s]->size / (_stacks[s]->recommendedChunksPerRank * _commSize));
				maximumWorkShare = (int(maximumWorkShare / _stacks[s]->recommendedChunkSizeMultiple) + 1) * _stacks[s]->recommendedChunkSizeMultiple;
				if (maximumWorkShare < 1) maximumWorkShare = 1;

				//determine dynamic chunk size according to compute power
				float myComputePower = _currentCalculationWorkDone[rank][s] / _currentCalculationTime[rank][s];
				float totalComputePower = 0;
				for (int i = 0; i < _commSize; ++i) 
					for (StackIdentifier j = 0; j < StackIdentifier(_stacks.size()); ++j) totalComputePower += _currentCalculationWorkDone[i][j] / _currentCalculationTime[i][j];
				
				int newCurrentCalculationProgress;
				if (std::isfinite(myComputePower) && std::isfinite(totalComputePower))
				{
					//factor in amount of remaining work
					int remainingWork = int(_stacks[s]->size - _currentCalculationStackProgress[s]);
					int myWorkShare = int((myComputePower / totalComputePower) * remainingWork);
					myWorkShare = (int(myWorkShare / _stacks[s]->recommendedChunkSizeMultiple) + 1) * _stacks[s]->recommendedChunkSizeMultiple;

					//clip to min/max chunk size
					int minumumWorkTime = 100;
					int minimumWorkShare = int(myComputePower * minumumWorkTime);
					if (minimumWorkShare < 1) minimumWorkShare = 1;
					if (myWorkShare < minimumWorkShare) myWorkShare = minimumWorkShare;
					if (myWorkShare > maximumWorkShare) myWorkShare = maximumWorkShare;
					newCurrentCalculationProgress = _currentCalculationStackProgress[s] + myWorkShare;
				}
				else newCurrentCalculationProgress = _currentCalculationStackProgress[s] + maximumWorkShare;

				//clip chunk to stack size
				if (newCurrentCalculationProgress >= _stacks[s]->size) newCurrentCalculationProgress = _stacks[s]->size;

				//write chunk parameters
				_currentCalculationWorkDone[rank][s] += newCurrentCalculationProgress - _currentCalculationStackProgress[s];
				c.properties[HMP_CHUNK_PROPERTY_STACK] = s;
				c.properties[HMP_CHUNK_PROPERTY_BEGIN] = decltype(c.properties[HMP_CHUNK_PROPERTY_BEGIN])(_currentCalculationStackProgress[s]);
				c.properties[HMP_CHUNK_PROPERTY_END] = decltype(c.properties[HMP_CHUNK_PROPERTY_END])(newCurrentCalculationProgress);
				Log::log << Log::LogLevel::Debug << "LoadManager spawned chunk (stack " << s << ", from " << _currentCalculationStackProgress[s] << ", to " << newCurrentCalculationProgress << ", rank  " << rank << ")" << Log::endl;

				//return chunk
				_currentCalculationStackProgress[s] = newCurrentCalculationProgress;
				break;
			}
			_currentCalculationChunkSpawntime[rank] = boost::posix_time::microsec_clock::local_time();
			_currentCalculationChunkSpawned[rank] = c;
			return c;
		}

		/**
		 * @brief Send a workload chunk to a specified MPI rank. 
		 * 
		 * @param rank Receiver's MPI rank. 
		 */
		void _issueChunk(const int rank)
		{
			#ifdef HMP_MPI_ENABLED
			Chunk c = _spawnChunk(rank);
			MPI_Send(&c.properties, 3, MPI_INT, rank, static_cast<int>(DataStackBase::MessageTag::Chunk), _communicator);

			if (!c.isVoid())
			{
				for (StackIdentifier i = 0; i < StackIdentifier(_stacks.size()); ++i)
				{
					//we may expect to receive data from additional slave stacks; due to MPI's non-overtaking policy, it is sufficient to only store the most recent request per rank
					if (i == c.properties[HMP_CHUNK_PROPERTY_STACK] || _stacks[i]->master == c.properties[HMP_CHUNK_PROPERTY_STACK]) _stacks[i]->receive(c.properties[HMP_CHUNK_PROPERTY_BEGIN], c.properties[HMP_CHUNK_PROPERTY_END] - c.properties[HMP_CHUNK_PROPERTY_BEGIN], rank, _communicator, _pendingRequests[rank]);
				}
			}
			else _pendingRequests[rank] = MPI_REQUEST_NULL;
			#endif
		}

		/**
		 * @brief Mark a workload chunk as completed. Relevant only for runtime statistics information. 
		 * 
		 * @param rank MPI rank which had completed the workload. 
		 */
		void _despawnChunk(const int rank)
		{
			float chunktime = float((boost::posix_time::microsec_clock::local_time() - _currentCalculationChunkSpawntime[rank]).total_milliseconds());
			std::lock_guard<std::mutex> lock(_currentCalculationChunkSpawnerLock);
			_currentCalculationTime[rank][_currentCalculationChunkSpawned[rank].properties[HMP_CHUNK_PROPERTY_STACK]] += chunktime;
		}

		float _totalCalculationTime; ///< Accumulated time in milliseconds which has been spent on calculate() calls over the lifetime of the LoadManager instance. 
		std::vector<float> *_totalComputeTime; ///< _totalComputeTime[rank][stack] is the accumulated time in milliseconds which MPI rank `rank` spent computing on `stack`. 
		std::vector<int> *_currentCalculationWorkDone; ///< _currentCalculationWorkDone[rank][stack] is the number of calculations which have been performed by MPI rank `rank` on `stack` in the current calculate() call. 
		std::vector<float> *_currentCalculationTime; ///< _currentCalculationTime[rank][stack] is the time in milliseconds spent by MPI rank `rank` until returning chunk result for `stack` in the current calculate() call. 
		boost::posix_time::ptime *_currentCalculationChunkSpawntime; ///< _currentCalculationChunkSpawntime[rank] specifies time at which the most recent chunk has been issued to MPI rank `rank` in the current calculate() call. 
		Chunk *_currentCalculationChunkSpawned; ///< _currentCalculationChunkSpawned[rank] stores the most recent chunk generated for MPI rank `rank` in the current calculate() call. 
		std::vector<float> _currentCalculationComputeTimeBuffer; ///< _currentCalculationComputeTimeBuffer[rank*_stacks.size()+stack] is a buffer for the time in milliseconds spent on computing `stack` in the current calculate() call. 
		std::vector<bool> _currentCalculationStackMask; ///< _currentCalculationStackMask[stack] specifies whether `stack` should be computed in the current calculate() call. 
		std::vector<int> _currentCalculationStackProgress; ///< _currentCalculationStackProgress[stack] specifies the current progress (pointer to the next unissued value) which has already been issued for computation in the current calculate() call. 
		std::mutex _currentCalculationChunkSpawnerLock; ///< Lock to synchronize chunk spawning for remote calculations and for local worker threads. 

		HMP_ENABLE_IF_MPI(MPI_Request *_pendingRequests); ///< MPI request objects associated with the return values for workload chunks that have been issued. 
	};

	/**
	 * @brief LoadManager slave implementation, responsible for receiving and executing workload chunks form a LoadManagerMaster implementation. 
	 * 
	 * @see LoadManager
	 */
	class LoadManagerSlave : public LoadManager
	{
		friend LoadManager *newLoadManager(const int serverRank HMP_APPEND_IF_MPI(const MPI_Comm communicator));
	public:
		/**
		 * @brief Calculate a list of stacks, where the stack identifiers are provided in list form. 
		 * 
		 * @param stackIds Pointer to the first StackIdentifier. 
		 * @param size Number of stacks. 
		 */
		virtual void calculate(const StackIdentifier *stackIds, const int size) override
		{
			#ifdef HMP_MPI_ENABLED
			memset(_currentCalculationComputeTimeBuffer.data(), 0, _currentCalculationComputeTimeBuffer.size() * sizeof(float));

			Chunk c;
			for (;;)
			{
				//receive chunk
				_waitChunk(c);
				if (c.isVoid()) break;

				//compute chunk
				boost::posix_time::ptime tic = boost::posix_time::microsec_clock::local_time();
				_calculateChunk(c);
				boost::posix_time::ptime toc = boost::posix_time::microsec_clock::local_time();
				_currentCalculationComputeTimeBuffer[c.properties[HMP_CHUNK_PROPERTY_STACK]] += (toc - tic).total_milliseconds();

				//return chunk
				_returnChunk(c);
			}

			//broadcast result
			for (int s = 0; s < size; ++s)
			{
				if (_stacks[stackIds[s]]->autoBroadcast) broadcast(stackIds[s]);
			}

			//gather compute time statistics
			MPI_Gather(_currentCalculationComputeTimeBuffer.data(), int(_currentCalculationComputeTimeBuffer.size()), MPI_FLOAT, nullptr, int(_currentCalculationComputeTimeBuffer.size()), MPI_FLOAT, _serverRank, _communicator);
			#endif
		}

	protected:
		/**
		 * @brief Construct a new LoadManagerSlave object
		 * 
		 * @param serverRank The MPI rank to take on the master role. 
		 * @param communicator The MPI communicator used for communication. 
		 */
		LoadManagerSlave(const int serverRank HMP_APPEND_IF_MPI(const MPI_Comm communicator)) : LoadManager(serverRank HMP_APPEND_IF_MPI(communicator)) {}

		/**
		 * @brief Register a DataStackBase with the LoadManager. By registering the stack, the LoadManager assumes responsibility to synchronize data between MPI ranks as necessary. 
		 * 
		 * @param stack The stack to be registered. 
		 * @return StackIdentifier Unique identifier of the newly registered stack; Used e.g. for LoadManager::calculate calls. 
		 */
		virtual StackIdentifier _registerStack(DataStackBase *stack) override
		{
			StackIdentifier identifier = LoadManager::_registerStack(stack);

			_currentCalculationComputeTimeBuffer.resize(_stacks.size());

			return identifier;
		}

		/**
		 * @brief Wait until a workload chunk has been received from a LoadManagerMaster instance. 
		 * 
		 * @param[out] chunk The workload definition which has been received. 
		 */
		void _waitChunk(Chunk &chunk) const
		{
			#ifdef HMP_MPI_ENABLED
			MPI_Status status;
			MPI_Recv(&chunk.properties, 3, MPI_INT, _serverRank, static_cast<int>(DataStackBase::MessageTag::Chunk), _communicator, &status);
			#endif
		}

		/**
		 * @brief Return the result of the workload defined by a specific chunk. 
		 * 
		 * @param chunk The workload definition whose results are to be returned. 
		 */
		void _returnChunk(const Chunk &chunk) const
		{
			#ifdef HMP_MPI_ENABLED
			for (int i = 0; i < int(_stacks.size()); ++i)
			{
				if (i == chunk.properties[HMP_CHUNK_PROPERTY_STACK] || _stacks[i]->master == chunk.properties[HMP_CHUNK_PROPERTY_STACK]) _stacks[i]->send(chunk.properties[HMP_CHUNK_PROPERTY_BEGIN], chunk.properties[HMP_CHUNK_PROPERTY_END] - chunk.properties[HMP_CHUNK_PROPERTY_BEGIN], _serverRank, _communicator);
			}
			#endif
		}

		std::vector<float> _currentCalculationComputeTimeBuffer; ///< _currentCalculationComputeTimeBuffer[stack] is a buffer for the time in milliseconds spent on computing `stack` in the current calculate() call. 
	};

	/**
	 * @brief Create a new LoadManager instance. 
	 * If the current MPI rank is the designated master rank, returns a LoadManagerMaster instance, otherwise returns a LoadManagerSlave instance. 
	 * 
	 * @param serverRank MPI rank to assume the master role. 
	 * @param communicator The MPI communicator to operate on. 
	 * @return LoadManager* Newly generated LoadManager instance. 
	 */
	inline LoadManager *newLoadManager(const int serverRank = 0 HMP_APPEND_IF_MPI(const MPI_Comm communicator = MPI_COMM_WORLD))
	{
		int myRank;
		HMP_ENABLE_IF_MPI(MPI_Comm_rank(communicator, &myRank));
		HMP_DISABLE_IF_MPI(myRank = 0);

		if (myRank == serverRank) return new LoadManagerMaster(serverRank HMP_APPEND_IF_MPI(communicator));
		else return new LoadManagerSlave(serverRank HMP_APPEND_IF_MPI(communicator));
	};

} //namespace HMP

#undef HMP_CHUNK_PROPERTY_STACK
#undef HMP_CHUNK_PROPERTY_BEGIN
#undef HMP_CHUNK_PROPERTY_END

#undef HMP_MPI_ENABLED
#undef HMP_ENABLE_IF_MPI
#undef HMP_DISABLE_IF_MPI
#undef HMP_APPEND_IF_MPI