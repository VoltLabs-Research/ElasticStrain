#pragma once

#include <opendxa/core/opendxa.h>
#include <opendxa/core/lammps_parser.h>
#include <nlohmann/json.hpp>
#include <opendxa/core/particle_property.h>
#include <opendxa/structures/crystal_structure_types.h>
#include <opendxa/analysis/structure_analysis.h>
#include <string>

namespace OpenDXA{
using json = nlohmann::json;
    
class ElasticStrainWrapper{
public:
    ElasticStrainWrapper();

    void setInputCrystalStructure(LatticeStructureType structure);
    void setIdentificationMode(StructureAnalysis::Mode mode);
    void setRMSD(float rmsd);

    void setParameters(
        double latticeConstant,
        double caRatio,
        bool pushForward,
        bool calculateDeformationGradient,
        bool calculateStrainTensors
    );

    json compute(
        const LammpsParser::Frame& frame,
        const std::string& outputFilename = ""
    );

private:
    LatticeStructureType _inputCrystalStructure;
    StructureAnalysis::Mode _identificationMode;
    float _rmsd;

    double _latticeConstant;
    double _caRatio;
    bool _pushForward;
    bool _calculateDeformationGradient;
    bool _calculateStrainTensors;


    std::shared_ptr<Particles::ParticleProperty> createPositionProperty(const LammpsParser::Frame &frame);
};

}