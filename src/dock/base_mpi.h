//
#ifndef BASE_MPI_H
#define BASE_MPI_H 

#include <string>
#include <vector>
class DOCKMol;
class Parameter_Reader;
class RANKMol;
#ifdef BUILD_DOCK_WITH_MPI
#  include "mpi.h"
#endif


class           Base_MPI {

  public:
    int             rank;       // machine rank
    int             comm_size;  // size of virtual machine
    int             mpi_work_unit;      // default 1

    std::vector < DOCKMol > send_queue;      // queue of mols to be sent via MPI
    std::vector < RANKMol > recv_queue;      // queue of mols received via MPI
    bool            receive_flag;       // used while looping over received
                                        // mols from clients

    bool            continue_mols;      // default true

    int             request[2]; // node rank; req type:0=no request, 1=send
                                // request, 2=recv request
    int             response[2];        // node rank; response type: 0=no
                                        // mols(shut down), n=# mols being sent
    int             stats[3];   // statistics on processed mols to be returned
                                // to server
    int             active_clients;     // number of active clients (only used
                                        // for master node)

    void            input_parameters(Parameter_Reader & parm);
    bool            initialize_mpi(int *argc, char ***argv);
    void            finalize_mpi();

    bool            continue_reading_mols();
    bool            is_master_node();
    bool            is_client_node();

    void            add_to_send_queue(DOCKMol &);
    bool            get_from_send_queue(DOCKMol &);

    void            add_to_recv_queue(RANKMol &);
    bool            get_from_recv_queue(RANKMol &);

    bool            listen_for_request();
    bool            is_receive_request();
    bool            is_send_request();

    bool            send_queue_request();
    void            transmit_send_queue();

    void            recv_queue_request();
    void            receive_recv_queue();

    void            barrier();

};

void            dockmol_to_string(DOCKMol &, std::string &);
void            string_to_dockmol(DOCKMol &, std::string &);

#endif  // BASE_MPI_H

