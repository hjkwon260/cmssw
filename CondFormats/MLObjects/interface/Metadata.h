#ifndef CondFormats_MLObjects_Metadata_h
#define CondFormats_MLObjects_Metadata_h

#include <string>
#include "CondFormats/Serialization/interface/Serializable.h"

class Metadata {
public:
    Metadata() {}
    Metadata(std::string model_name, int version) : model_name_(model_name), version_(version) {}

    std::string model_name() const { return model_name_; }
    int version() const { return version_; }

private:

    std::string model_name_;
    int version_;

    COND_SERIALIZABLE;
};

class MetadataCollection {
public:
    MetadataCollection() {}

    void add_model(const Metadata& m) { models_.push_back(m); }
    const std::vector<Metadata>& models() const { return models_; }

private:
    std::vector<Metadata> models_;

    COND_SERIALIZABLE;
};

#endif
