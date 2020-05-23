// Author: Armin Töpfer

#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <pbbam/FastaReader.h>
#include <pbbam/FastaWriter.h>
#include <pbbam/FastqReader.h>
#include <pbbam/FastqWriter.h>

#include <pbcopper/third-party/edlib.h>
#include <pbcopper/utility/SequenceUtils.h>

namespace {
struct EdlibAlignResultHandle
{
    EdlibAlignResultHandle(EdlibAlignResult aln) : data_(std::move(aln)) {}

    ~EdlibAlignResultHandle() { edlibFreeAlignResult(data_); }

    EdlibAlignResult data_;
};

EdlibAlignConfig DefaultEdlibAlignConfig()
{
    static const EdlibAlignConfig config =
        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, nullptr, 0);
    return config;
}

EdlibAlignResultHandle EdlibAlign(const char* target, const int targetLength, const char* query,
                                  const int queryLength, const EdlibAlignConfig& config)
{
    // edlib flips the params here
    return EdlibAlignResultHandle{edlibAlign(query, queryLength, target, targetLength, config)};
}

bool AlignsEdlib(const std::string& refSeq, const std::string& qry, const int32_t qrySize,
                 const double minIdentityPerc)
{
    EdlibAlignConfig config = EdlibAlignConfig(DefaultEdlibAlignConfig());
    // Inform edlib of the maximum edit distance that is interesting.
    // config.k = 1 + (100.0 - minIdentityPerc) * (qrySize / 100.0);

    const auto aln = EdlibAlign(refSeq.c_str(), refSeq.size(), qry.c_str(), qrySize, config);
    // Check edit distance (-1 if Edlib does not find an alignment)
    if (aln.data_.editDistance == -1) return false;

    const double identityPerc =
        100.0 * (1.0 - (1.0 * aln.data_.editDistance / aln.data_.alignmentLength));
    const bool pass = identityPerc >= minIdentityPerc;
    std::cerr << identityPerc << ' ' << aln.data_.alignmentLength << ' ' << aln.data_.editDistance
              << ' ' << pass << '\n';
    return identityPerc >= minIdentityPerc;
}
}  // namespace

template <typename R, typename W>
void Process(const std::string& input, const std::string& outputFiltered,
             const std::string& outputLacking, const std::vector<std::string>& primers,
             const double minIdentity)
{
    W fastqFiltered{outputFiltered};
    W fastqLacking{outputLacking};
    int32_t filtered = 0;
    int32_t lacking = 0;
    for (const auto& record : R{input}) {
        bool found = false;
        for (const auto& p : primers) {
            found = AlignsEdlib(record.Bases(), p, p.size(), minIdentity);
            if (found) break;
        }

        if (found) {
            ++filtered;
            fastqFiltered.Write(record);
        } else {
            ++lacking;
            fastqLacking.Write(record);
        }
    }
    std::cout << "FILTERED : " << filtered << " (" << (100.0 * filtered / (filtered + lacking))
              << "%)\n"
              << "NO HITS  : " << lacking << '\n';
}

int main(int argc, char* argv[])
{
    using namespace PacBio;

    const std::string input{argv[1]};
    const std::string outputFiltered{argv[2]};
    const std::string outputLacking{argv[3]};
    const std::string primer{argv[4]};
    int32_t minIdentity = 90;
    if (argc == 6) {
        minIdentity = std::stoi(std::string{argv[5]});
        std::cerr << "Overriding default minimum identity to " << minIdentity << '\n';
    }

    std::vector<std::string> primers;
    for (const auto& record : BAM::FastaReader{primer}) {
        primers.emplace_back(record.Bases());
        primers.emplace_back(Utility::ReverseComplemented(record.Bases()));
    }
    if (boost::iends_with(input, ".fa") || boost::iends_with(input, ".fasta") ||
        boost::iends_with(input, ".fa.gz") || boost::iends_with(input, ".fasta.gz")) {
        Process<BAM::FastaReader, BAM::FastaWriter>(input, outputFiltered, outputLacking, primers,
                                                    minIdentity);
    } else if (boost::iends_with(input, ".fq") || boost::iends_with(input, ".fastq") ||
               boost::iends_with(input, ".fq.gz") || boost::iends_with(input, ".fastq.gz")) {
        Process<BAM::FastqReader, BAM::FastqWriter>(input, outputFiltered, outputLacking, primers,
                                                    minIdentity);
    }
}
