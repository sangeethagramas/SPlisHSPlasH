#include "Elasticity_Custom_Boundary.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include <extern/md5/md5.h>
#include "Utilities/BinaryFileReaderWriter.h"
#include "Utilities/StringTools.h"
#include <Utilities/FileSystem.h>
#include <cnpy.h>
#include <cstdint>
#include <array>

using namespace SPH;
using namespace GenParam;

int SPH::Elasticity_Custom_Boundary::BOUNDARY_DATA_FILE = -1;

Elasticity_Custom_Boundary::Elasticity_Custom_Boundary(FluidModel *model) :
	Elasticity_Kugelstadt2021(model)
{
	std::cout << "Elasticity_Custom_Boundary::Elasticity_Custom_Boundary" << std::endl;
}

Elasticity_Custom_Boundary::~Elasticity_Custom_Boundary(void)
{
	std::cout << "Elasticity_Custom_Boundary::~Elasticity_Custom_Boundary" << std::endl;
}

void Elasticity_Custom_Boundary::initParameters()
{
	Elasticity_Kugelstadt2021::initParameters();
	std::cout << "Elasticity_Custom_Boundary::initParameters" << std::endl;

    ParameterBase::GetFunc<std::string> getBdyDF = [&]()-> std::string { return m_boundary_data_file; };
	ParameterBase::SetFunc<std::string> setBdyDF = [&](std::string val) {
        m_boundary_data_file = val;
    };
    BOUNDARY_DATA_FILE = createStringParameter("boundaryDataFile", "boundary data filepath", getBdyDF, setBdyDF);
    setGroup(BOUNDARY_DATA_FILE, "Fluid Model|Elasticity");
	setDescription(BOUNDARY_DATA_FILE, "Path to .npz file with data");
	getParameter(BOUNDARY_DATA_FILE)->setReadOnly(true);

}


void Elasticity_Custom_Boundary::determineFixedParticles()
{
	// If no file is provided, keep all particles active.
	if (m_boundary_data_file.empty())
	{
		LOG_WARN << "Elasticity_Custom_Boundary: No boundary data file set; skipping fixed particle assignment.";
		return;
	}

	std::cout << "Elasticity_Custom_Boundary::determineFixedParticles: m_boundary_data_file = " << m_boundary_data_file << std::endl;

	// Ensure the file exists before attempting to load.
	if (!Utilities::FileSystem::fileExists(m_boundary_data_file))
	{
		LOG_ERR << "Elasticity_Custom_Boundary: Boundary data file not found: " << m_boundary_data_file;
		return;
	}

	try
	{
		// reset recorded order of fixed boundary indices
		m_bdy_indices.clear();
		m_bdy_num_boundary = 0u;

		// Load only the 'boundary_indices' entry from the NPZ
		cnpy::NpyArray arr = cnpy::npz_load(m_boundary_data_file, "boundary_indices");

		// Determine number of elements (product of shape)
		size_t count = 1;
		for (size_t d : arr.shape)
			count *= d;
		if (count == 0)
		{
			LOG_WARN << "Elasticity_Custom_Boundary: 'boundary_indices' is empty in " << m_boundary_data_file;
			return;
		}

		const unsigned int numParticles = m_model->numActiveParticles();
		size_t fixedCount = 0;

		auto mark_if_valid = [&](size_t idx)
		{
			if (idx < numParticles)
			{
				m_model->setParticleState(static_cast<unsigned int>(idx), ParticleState::Fixed);
				m_bdy_indices.push_back(static_cast<unsigned int>(idx));		// preserve order of marking
				fixedCount++;
			}
			else
			{
				LOG_WARN << "Elasticity_Custom_Boundary: Out-of-range boundary index " << idx << " (numParticles=" << numParticles << ")";
			}
		};

		// Support common integer dtypes by checking word_size
		if (arr.word_size == sizeof(int64_t))
		{
			const int64_t* data = arr.data<int64_t>();
			for (size_t i = 0; i < count; i++)
				if (data[i] >= 0) mark_if_valid(static_cast<size_t>(data[i]));
		}
		else if (arr.word_size == sizeof(uint64_t))
		{
			const uint64_t* data = arr.data<uint64_t>();
			for (size_t i = 0; i < count; i++)
				mark_if_valid(static_cast<size_t>(data[i]));
		}
		else if (arr.word_size == sizeof(int32_t))
		{
			const int32_t* data = arr.data<int32_t>();
			for (size_t i = 0; i < count; i++)
				if (data[i] >= 0) mark_if_valid(static_cast<size_t>(data[i]));
		}
		else if (arr.word_size == sizeof(uint32_t))
		{
			const uint32_t* data = arr.data<uint32_t>();
			for (size_t i = 0; i < count; i++)
				mark_if_valid(static_cast<size_t>(data[i]));
		}
		else
		{
			LOG_ERR << "Elasticity_Custom_Boundary: Unsupported dtype (word_size=" << arr.word_size << ") for 'boundary_indices' in " << m_boundary_data_file;
			return;
		}

		LOG_INFO << "Elasticity_Custom_Boundary: Marked " << fixedCount << " particles as Fixed from 'boundary_indices'.";
		m_bdy_num_boundary = static_cast<unsigned int>(m_bdy_indices.size());

		// Also load the full boundary position sequence once and cache it
		{
			cnpy::NpyArray arr_seq = cnpy::npz_load(m_boundary_data_file, "boundary_pos_seq");
			if (arr_seq.shape.size() != 3)
			{
				LOG_ERR << "Elasticity_Custom_Boundary: 'boundary_pos_seq' must be 3D array [num_frames, num_boundary, 3].";
				return;
			}
			const size_t F = arr_seq.shape[0];
			const size_t B = arr_seq.shape[1];
			const size_t C = arr_seq.shape[2];
			if (C != 3)
			{
				LOG_ERR << "Elasticity_Custom_Boundary: Last dim of 'boundary_pos_seq' must be 3.";
				return;
			}
			if (m_bdy_num_boundary != 0u && B != m_bdy_num_boundary)
			{
				LOG_WARN << "Elasticity_Custom_Boundary: 'boundary_pos_seq' boundary count (" << B << ") does not match number of fixed boundary particles (" << m_bdy_num_boundary << ").";
			}

			m_bdy_num_frames = static_cast<unsigned int>(F);
			if (m_bdy_num_boundary == 0u)
				m_bdy_num_boundary = static_cast<unsigned int>(B);
			m_bdy_pos_seq.resize(F * B);

			// Copy and convert to Vector3r
			if (arr_seq.word_size == sizeof(double))
			{
				const double* d = arr_seq.data<double>();
				for (size_t f = 0; f < F; ++f)
				{
					for (size_t b = 0; b < B; ++b)
					{
						const size_t base = (f * B + b) * 3;
						Vector3r p;
						p[0] = static_cast<Real>(d[base + 0]);
						p[1] = static_cast<Real>(d[base + 1]);
						p[2] = static_cast<Real>(d[base + 2]);
						m_bdy_pos_seq[f * B + b] = p;
					}
				}
			}
			else if (arr_seq.word_size == sizeof(float))
			{
				const float* d = arr_seq.data<float>();
				for (size_t f = 0; f < F; ++f)
				{
					for (size_t b = 0; b < B; ++b)
					{
						const size_t base = (f * B + b) * 3;
						Vector3r p;
						p[0] = static_cast<Real>(d[base + 0]);
						p[1] = static_cast<Real>(d[base + 1]);
						p[2] = static_cast<Real>(d[base + 2]);
						m_bdy_pos_seq[f * B + b] = p;
					}
				}
			}
			else
			{
				LOG_ERR << "Elasticity_Custom_Boundary: Unsupported dtype for 'boundary_pos_seq'.";
				return;
			}
			m_bdy_loaded = true;
			m_bdy_frame_index = 0u;
			LOG_INFO << "Elasticity_Custom_Boundary: Cached boundary animation (frames=" << m_bdy_num_frames << ", boundary=" << m_bdy_num_boundary << ").";
		}
	}
	catch (const std::exception& e)
	{
		LOG_ERR << "Elasticity_Custom_Boundary: Failed to load '" << m_boundary_data_file << "': " << e.what();
	}
}

void Elasticity_Custom_Boundary::updateBoundaryParticles()
{
	// Use cached data only
	if (!m_bdy_loaded || m_bdy_num_frames == 0 || m_bdy_num_boundary == 0)
		return;
	
	if (m_bdy_frame_index >= m_bdy_num_frames)
		return;

	// Select frame (loop)
	const unsigned int f = m_bdy_frame_index;

	const unsigned int numParticles = m_model->numActiveParticles();
	const size_t B = static_cast<size_t>(m_bdy_num_boundary);
	for (size_t k = 0; k < B; ++k)
	{
		const unsigned int pid = m_bdy_indices[k];
		if (pid >= numParticles) 
			continue;
		if (m_model->getParticleState(pid) == ParticleState::Fixed)
		{
			Vector3r& xi = m_model->getPosition(pid);
			xi = m_bdy_pos_seq[f * B + k];
		}
	}

	// Advance frame counter
	m_bdy_frame_index++;
}


void Elasticity_Custom_Boundary::step()
{
    updateBoundaryParticles();
	// apply accelerations
	Elasticity_Kugelstadt2021::step();
}



