#ifndef CondFormats_MLObjects_MetadataWrapper_h
#define CondFormats_MLObjects_MetadataWrapper_h

#include <string>
#include "CondFormats/Serialization/interface/Serializable.h"

class MetadataWrapper {
public:
    MetadataWrapper() = default;
    MetadataWrapper(std::string model_name, std::string version,
                    std::string model_path, std::string preproc_path)
        : model_name_(model_name), version_(version),
          model_path_(model_path), preproc_path_(preproc_path) {}

    std::string model_name() const { return model_name_; }
    std::string version() const { return version_; }
    std::string model_path() const { return model_path_; }
    std::string preprocessing_path() const { return preproc_path_; }

private:
    std::string model_name_;
    std::string version_;
    std::string model_path_;
    std::string preproc_path_;

    COND_SERIALIZABLE;

};

class MetadataWrapperCollection {
public:
    MetadataWrapperCollection() = default;

    void add(const MetadataWrapper& m) { models_.push_back(m); }
    const std::vector<MetadataWrapper>& models() const { return models_; }

private:
    std::vector<MetadataWrapper> models_;
    COND_SERIALIZABLE;

};

#endif

