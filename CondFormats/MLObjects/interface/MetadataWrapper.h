#ifndef CondFormats_MLObjects_MetadataWrapper_h
#define CondFormats_MLObjects_MetadataWrapper_h

#include <string>
#include "CondFormats/Serialization/interface/Serializable.h"

class MetadataWrapper {
public:
    MetadataWrapper() = default;
    MetadataWrapper(std::string name, std::string version,
                    std::string modelPath, std::string preprocPath)
        : name_(name), version_(version),
          modelPath_(modelPath), preprocPath_(preprocPath) {}

    std::string name() const { return name_; }
    std::string version() const { return version_; }
    std::string modelPath() const { return modelPath_; }
    std::string preprocessingPath() const { return preprocPath_; }

private:
    std::string name_;
    std::string version_;
    std::string modelPath_;
    std::string preprocPath_;

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

