#ifndef CondFormats_MLObjects_Metadata_h
#define CondFormats_MLObjects_Metadata_h

#include <string>
#include "CondFormats/Serialization/interface/Serializable.h"

class Metadata {
public:
    Metadata() : value_(0) {}
    Metadata(int v, std::string t) : value_(v), info_(t) {}

    int value() const { return value_; }
    std::string info() const { return info_; }

private:
    int value_;
    std::string info_;

    COND_SERIALIZABLE;
};

class MetadataCollection {
public:
    MetadataCollection() {}

    void addModel(const Metadata& m) { models_.push_back(m); }
    const std::vector<Metadata>& models() const { return models_; }

private:
    std::vector<Metadata> models_;

    COND_SERIALIZABLE;
};

#endif
