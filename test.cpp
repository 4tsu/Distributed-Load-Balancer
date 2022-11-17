void export_cdview(std::vector<int> atoms, int rank) {
    static int count = 0;
    char filename[256];
    sprintf(filename, "conf%03d.cdv", count);
    ++count;
    std::ofstream ofs(filename, std::ios::app);
    if (rank==0) {
        ofs << "Header" << std::endl;
        ofs << "Header" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
		for (auto a : atoms) {
        ofs << rank << " ";
        ofs << a << " ";
        ofs << std::endl;
		}
}



void export_cdview_independent(std::vector<int> atoms, int rank) {
    static int count = 0;
    char filename[256];
    sprintf(filename, "tconf%03d_%d.cdvt", count, rank);
    ++count;
    std::ofstream ofs(filename, std::ios::app);
    if (rank==0) {
        ofs << "Header" << std::endl;
        ofs << "Header" << std::endl;
    }
		for (auto a : atoms) {
        ofs << rank << " ";
        ofs << a << " ";
        ofs << std::endl;
		}
}



void concatenate_cdview(int rank, int procs) {
		int count = 0;
		char output[256];
		sprintf(output, "conf%03d.cdv", count);
    std::ofstream ofs(output, std::ios::app);
		for (int i=0; i<procs; i++) {
				char filename[256];
				sprintf(filename, "tconf%03d_%d.cdvt", count, i);
				std::ifstream reading_file;
				reading_file.open(filename, std::ios::in);
				std::string line;
				while(std::getline(reading_file, line)) {
					ofs << line << std::endl;
				}
		}
#ifdef FS
    if (mi.rank == 0) {
        for (const auto & file : std::filesystem::directory_iterator("./")) {
            std::string path = file.path();
            size_t word_pos = path.find(".cdvt");
            if (word_pos != std::string::npos) {
                std::filesystem::remove(path);
                continue;
            }
        }
    }
#endif
}



int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	int rank, procs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	std::vector<int> v(10000);
	std::fill(v.begin(), v.end(), rank);
	// export_cdview(v, rank);
	export_cdview_independent(v, rank);
	if (rank==0) concatenate_cdview(rank, procs);
	MPI_Finalize();
	return 0;
}
