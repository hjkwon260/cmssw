/**
 * Author: H. Kwon
 */

#ifndef CondFormats_MLObjects_MetadataWrapper_h
#define CondFormats_MLObjects_MetadataWrapper_h

#include <string>
#include <vector>
#include <memory>
#include "CondFormats/Serialization/interface/Serializable.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

class MetadataWrapper {
public:
    MetadataWrapper() = default;
    MetadataWrapper(std::string model_name, std::string version,
                    std::string hash, std::string model_path, std::string preproc_path)
        : model_name_(model_name), version_(version), hash_(hash),
          model_path_(model_path), preproc_path_(preproc_path) {}

    std::string model_name() const { return model_name_; }
    std::string version() const { return version_; }
    std::string hash() const { return hash_; }
    std::string model_path() const { return model_path_; }
    std::string preprocessing_path() const { return preproc_path_; }
    const std::vector<std::string>& input_features() const { return input_features_; }
    const std::vector<std::string>& flav_names() const { return flav_names_; }

    void setOnnxRuntime(std::shared_ptr<cms::Ort::ONNXRuntime> runtime) {
        runtime_ = runtime;
    }

    std::shared_ptr<cms::Ort::ONNXRuntime> onnxRuntime() const { return runtime_; }

    void set_input_features(const std::vector<std::string>& features) { input_features_ = features; }
    void set_flav_names(const std::vector<std::string>& names) { flav_names_ = names; }


private:
    std::string model_name_;
    std::string version_;
    std::string hash_;
    std::string model_path_;
    std::string preproc_path_;
    std::vector<std::string> input_features_;
    std::vector<std::string> flav_names_;
    std::shared_ptr<cms::Ort::ONNXRuntime> runtime_;

    // COND_SERIALIZABLE;

};

#endif

