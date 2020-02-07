// Author: Armin Töpfer

#include <algorithm>
#include <cstdlib>
#include <random>
#include <string>

#include <pbbam/BamHeader.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>
#include <pbbam/Compare.h>
#include <pbbam/PbiRawData.h>
#include <pbbam/ZmwGroupQuery.h>
#include <pbcopper/data/LocalContextFlags.h>

std::vector<int32_t> UniqueZmws(const PacBio::BAM::DataSet ds, bool dieOnError = true)
{
    using namespace PacBio::BAM;
    const auto bamFiles = ds.BamFiles();
    if (bamFiles.size() != 1) {
        std::cerr << "Chunking only works with one input BAM file!" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (!bamFiles[0].PacBioIndexExists()) {
        std::cerr << "PBI file is missing for input BAM file! Please create one using pbindex!"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    bool fromZMWSet = false;
    bool toZMWSet = false;
    int32_t fromZMW{};
    int32_t toZMW{};
    using Operator = PacBio::BAM::Compare::Type;
    for (const auto& filter : ds.Filters()) {
        for (const auto& property : filter.Properties()) {
            int32_t value = -1;
            Operator op = Operator::NOT_CONTAINS;
            bool isZMW = false;
            for (const auto& attribute : property.Attributes()) {
                if (attribute.first == "Name" && attribute.second == "zm") {
                    isZMW = true;
                } else if (attribute.first == "Operator" && attribute.second == "<") {
                    op = Operator::LESS_THAN;
                } else if (attribute.first == "Operator" && attribute.second == "<=") {
                    op = Operator::LESS_THAN_EQUAL;
                } else if (attribute.first == "Operator" && attribute.second == ">") {
                    op = Operator::GREATER_THAN;
                } else if (attribute.first == "Operator" && attribute.second == ">=") {
                    op = Operator::GREATER_THAN_EQUAL;
                } else if (attribute.first == "Value") {
                    bool isNumber = true;
                    for (const auto c : attribute.second) {
                        if (!std::isdigit(c)) {
                            isNumber = false;
                            break;
                        }
                    }
                    if (isNumber) value = std::stoi(attribute.second);
                }
            }

            if (!isZMW) continue;

            switch (op) {
                case Operator::LESS_THAN_EQUAL:
                    toZMW = value + 1;
                    toZMWSet = true;
                    break;
                case Operator::LESS_THAN:
                    toZMW = value;
                    toZMWSet = true;
                    break;
                case Operator::GREATER_THAN:
                    fromZMW = value;
                    fromZMWSet = true;
                    break;
                case Operator::GREATER_THAN_EQUAL:
                    fromZMW = value - 1;
                    fromZMWSet = true;
                    break;
                default:
                    std::cerr << "WARN: Unset operator type for XML filter\n";
                    break;
            }
        }
    }
    const PbiRawData index{bamFiles[0].PacBioIndexFilename()};
    const auto& basicData = index.BasicData();
    const auto& zmws = basicData.holeNumber_;
    std::vector<int32_t> zmwsUniq;
    const int32_t numRecords = zmws.size();
    zmwsUniq.reserve(numRecords);
    if (numRecords == 0) {
        std::cerr << "No input records in PBI file!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (fromZMWSet && toZMWSet) {
        std::cerr << "ZMW filter range (" << fromZMW << ',' << toZMW << ')' << std::endl;
        if (zmws[0] > fromZMW && zmws[0] < toZMW) zmwsUniq.emplace_back(zmws[0]);
        for (int32_t i = 1; i < numRecords; ++i) {
            if (zmws[i] != zmws[i - 1] && zmws[i] > fromZMW && zmws[i] < toZMW)
                zmwsUniq.emplace_back(zmws[i]);
        }
    } else {
        zmwsUniq.emplace_back(zmws[0]);
        for (int32_t i = 1; i < numRecords; ++i) {
            if (zmws[i] != zmws[i - 1]) zmwsUniq.emplace_back(zmws[i]);
        }
    }
    std::cerr << "UNIQUE ZMWS: " << zmwsUniq.size() << std::endl;
    return zmwsUniq;
}

int main(int argc, char* argv[])
{
    const std::string input{argv[1]};
    const std::string output{argv[2]};
    const int32_t percentage = std::atoi(argv[3]);
    using namespace PacBio::BAM;
    const auto GetHeader = [&]() {
        BamReader reader{input};
        return reader.Header().DeepCopy();
    };
    BamHeader header = GetHeader();

    std::vector<int32_t> uniques = UniqueZmws(input);
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(uniques.begin(), uniques.end(), g);

    int32_t numberOfZmws = std::round(uniques.size() * percentage / 100.0);
    std::vector<int32_t> subset(uniques.begin(), uniques.begin() + numberOfZmws);

    ZmwGroupQuery qry{subset, input};
    BamWriter writer{output, header};
    int32_t zmwCounter = 0;
    int32_t subreadCounter = 0;
    for (const auto& zmw : qry) {
        for (auto& record : zmw) {
            writer.Write(record);
            ++subreadCounter;
        }
        ++zmwCounter;
    }
    std::cerr << "Written ZMWs: " << zmwCounter << std::endl;
    std::cerr << "Written subreads: " << subreadCounter << std::endl;
}
