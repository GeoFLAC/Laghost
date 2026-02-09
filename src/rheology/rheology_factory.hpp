// Copyright (c) 2024, Laghost Authors
// SPDX-License-Identifier: BSD-3-Clause
//
// Rheology Factory
// Factory pattern with auto-registration for constitutive models

#ifndef LAGHOST_RHEOLOGY_FACTORY_HPP
#define LAGHOST_RHEOLOGY_FACTORY_HPP

#include "constitutive_model.hpp"
#include <map>
#include <memory>
#include <functional>
#include <string>
#include <stdexcept>

namespace mfem
{
namespace geodynamics
{

/// Factory for creating constitutive model instances
class RheologyFactory
{
public:
   using Creator = std::function<std::unique_ptr<ConstitutiveModel>()>;
   
   /// Register a model creator with a name
   static void Register(const std::string &name, Creator creator)
   {
      Registry()[name] = std::move(creator);
   }
   
   /// Create a model instance by name
   static std::unique_ptr<ConstitutiveModel> Create(const std::string &name)
   {
      auto it = Registry().find(name);
      if (it == Registry().end())
      {
         throw std::runtime_error("Unknown constitutive model: " + name + 
                                  ". Available: " + AvailableModels());
      }
      return it->second();
   }
   
   /// Check if a model is registered
   static bool IsRegistered(const std::string &name)
   {
      return Registry().find(name) != Registry().end();
   }
   
   /// Get list of available models
   static std::string AvailableModels()
   {
      std::string result;
      for (const auto &pair : Registry())
      {
         if (!result.empty()) result += ", ";
         result += pair.first;
      }
      return result;
   }

private:
   /// Singleton registry map
   static std::map<std::string, Creator>& Registry()
   {
      static std::map<std::string, Creator> registry;
      return registry;
   }
};

/// Auto-registration helper macro
/// Usage: REGISTER_RHEOLOGY(MohrCoulombModel, "mohr_coulomb")
#define REGISTER_RHEOLOGY(ClassName, Name) \
   namespace { \
      static bool _##ClassName##_registered = []() { \
         ::mfem::geodynamics::RheologyFactory::Register(Name, []() { \
            return std::make_unique<ClassName>(); \
         }); \
         return true; \
      }(); \
   }

} // namespace geodynamics
} // namespace mfem

#endif // LAGHOST_RHEOLOGY_FACTORY_HPP
