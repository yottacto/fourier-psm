#pragma once

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include "lib/cpptoml/cpptoml.hh"

namespace icesp
{
    namespace configuration
    {
        struct case_
        {
            using toml_type = std::shared_ptr<cpptoml::table>;

            case_(toml_type toml)
            {
                #define GET(TYPE, NAME) do { \
                        auto NAME ## _opt = toml->template get_as<TYPE>(#NAME); \
                        if (NAME ## _opt) { \
                            NAME = std::move(* NAME ## _opt); \
                        } else { \
                            throw std::runtime_error{#NAME " must be " #TYPE}; \
                        } \
                    } while(0)
                GET(int, bncd);
                GET(int, bnpcd);
                GET(int, iteration);
                #undef GET
            }

            int bncd;
            int bnpcd;
            int iteration;
        };

        struct config
        {
            using toml_type = std::shared_ptr<cpptoml::table>;
            using policy_type = cpptoml::option<std::string>;
            using selection_type = cpptoml::option<std::vector<std::string>>;

            config(toml_type toml)
            {
                auto policy = toml->get_qualified_as<std::string>("mask.policy");
                auto selection = toml->get_qualified_array_of<std::string>("mask.selection");
                auto filter = filter_from_policy(policy, selection);

                for (auto & case_ : *toml->get_table("case")) {
                    if (filter(case_.first))
                        cases.emplace(case_.first, case_.second->as_table());
                }
            }

            std::unordered_map<std::string, case_> cases;

        private:
            static auto filter_from_policy(
                    policy_type const& policy,
                    selection_type const& selection)
                -> std::function<bool (std::string const&)>
            {
                if (!policy) return [](auto&) { return true; };
                if (!selection) {
                    throw std::runtime_error{
                        "selection must be used with policy"
                    };
                }

                std::unordered_set<std::string> set{selection->begin(), selection->end()};
                if (*policy == "allow") {
                    return [set=std::move(set)](auto& name) {
                        return set.count(name) != 0;
                    };
                }
                if (*policy == "deny") {
                    return [set=std::move(set)](auto& name) {
                        return set.count(name) == 0;
                    };
                }

                throw std::runtime_error{
                    "invalid policy, must be one of: \"allow\", \"deny\""
                };
            }
        };
    }
}

