//
// Created by Katie Barr (EI) on 23/11/2017.
//

#include "PhaseScaffolder.h"



PhaseScaffolder::PhaseScaffolder(std::string gfa_filename): sg(SequenceGraph()), mapper(PairedReadMapper(sg)){
    sg.load_from_gfa(gfa_filename);

}

void PhaseScaffolder::load_mappings(std::string r1_filename, std::string r2_filename){
    mapper.map_reads(r1_filename, r2_filename, prm10x);
    std::cout << "Mapped " << mapper.read_to_node.size() << " reads to " <<  mapper.reads_in_node.size() << "nodes" << std::endl;
}

void PhaseScaffolder::phase_components() {

//find and phase each component of gfa
    auto components = sg.connected_components();
// this finds 2 components for test graph...
    std::cout << "Found " << components.size() << " connected components " << std::endl;
    for (auto component:components) {
        HaplotypeScorer hs(sg);
        for (auto n:component) {
            std::cout << sg.oldnames[n] << std::endl; //" ";
        }
        std::cout << std::endl;
// should
        auto bubbles = sg.find_bubbles(component);
        hs.find_possible_haplotypes(bubbles);

        //hs.count_barcode_votes(reads1, reads2);
        hs.decide_barcode_haplotype_support();
// now have mappings and barcode support
        if (hs.barcode_haplotype_mappings.size() > 0) {
            hs.score_haplotypes();
        }
    }
}