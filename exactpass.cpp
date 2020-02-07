// Author: Armin Töpfer

#include <cstdlib>
#include <string>

#include <pbbam/BamHeader.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>
#include <pbbam/ZmwGroupQuery.h>
#include <pbcopper/data/LocalContextFlags.h>

int main(int argc, char* argv[])
{
    const std::string input{argv[1]};
    const std::string output{argv[2]};
    const int32_t maxPasses = std::atoi(argv[3]);
    using namespace PacBio::BAM;
    const auto GetHeader = [&]() {
        BamReader reader{input};
        return reader.Header().DeepCopy();
    };
    BamHeader header = GetHeader();

    ZmwGroupQuery qry{input};
    BamWriter writer{output, header};
    for (const auto& zmw : qry) {
        int32_t cov = 0;
        std::vector<BamRecord> records;
        for (auto& record : zmw) {
            const auto cx = record.LocalContextFlags();
            if (cx & LocalContextFlags::ADAPTER_BEFORE && cx & LocalContextFlags::ADAPTER_AFTER) {
                ++cov;
                records.emplace_back(std::move(record));
            }
            if (cov == maxPasses) break;
        }
        if (static_cast<int32_t>(records.size()) == maxPasses) {
            for (const auto& r : records)
                writer.Write(r);
        }
    }
}
