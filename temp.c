csc *** Abl_t = (csc***)malloc( 2 * sizeof(csc**));
    Abl_t[0] = Abl;
    Abl_t[1] = (csc**)malloc( total_procs * sizeof(csc*) );
    csc **no_name = Abl_t[1];
    for (int i = 0; i < total_procs; i++)
        no_name[i] = initCsc(Abl[0]->rowS, Abl[0]->colS, Abl[0]->nnz);
    

    for (int i = 1; i < total_procs; i++)
    {
        // Send - (Receive) data to the next (from the previous)
        csc ** row_used = Abl_t[ (i-1)%2 ];
        csc ** row_transfered = Abl_t[ i%2 ];

        if(!proc_num) printf("used:%d \n", proc_num, i);

        for(int j=0; j < total_procs; j++){
            csc *transfer_block = row_transfered[j];
            csc *use_block = row_used[j];

            int buf[2] = {transfer_block->rowS, transfer_block->nnz};
            MPI_Isend(buf, 2, MPI_INT, RIGHT, TAG_PARAM, MPI_COMM_WORLD , &reqs[ j * 3 ]);
            MPI_Isend(transfer_block->col_ptr, transfer_block->rowS + 1, MPI_INT, RIGHT, TAG_COL, MPI_COMM_WORLD , &reqs[ j * 3 ] + 1 );
            MPI_Isend(transfer_block->row, transfer_block->nnz, MPI_INT, RIGHT, TAG_ROW, MPI_COMM_WORLD , &reqs[ j * 3 ] + 2 );

            printf("P%d Started sending block\n", proc_num, i);

            int buf_used[2];
            MPI_Recv(buf_used, 2, MPI_INT, LEFT, TAG_PARAM, MPI_COMM_WORLD, MPI_STATUS_IGNORE); use_block->nnz = buf_used[1]; use_block->nnz = buf_used[1];
            MPI_Recv(use_block->col_ptr, use_block->rowS + 1, MPI_INT, LEFT, TAG_COL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(use_block->row, use_block->nnz, MPI_INT, LEFT, TAG_ROW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            printf("P%d Received block\n", proc_num, i);
        }



        // block_bmm(row_used, Bbl, Cbl[i], temps,total_procs,proc_num);

        // MPI_Waitall(3 * total_procs,reqs, MPI_STATUSES_IGNORE);

        printf("P%d i=%d\n", proc_num, i);
    }









/*for (int i = 1; i < total_procs; i++)
    {

        //   MPI_Waitall(3 * total_procs,reqs, MPI_STATUSES_IGNORE);

        // Send data to the next
        for (int s = 0; s < total_procs; s++)
            isend_block(Bbl[s], total_procs, RIGHT, &reqs[3 * s]);

        Bbl = rcv_blocks(total_procs, LEFT);
        printf("P%d Received Bbl\n", proc_num);

        for (int s = 0; s < total_procs; s++)
            isend_block(Abl[s], total_procs, RIGHT, &reqs[3 * s]);

        Abl = rcv_blocks(total_procs, LEFT);
        printf("P%d Received Abl\n", proc_num);

        //   block_bmm(Abl, Bbl, Cbl[i], temps,total_procs,proc_num);

        if (!proc_num)
            printCsc(Abl[0]);

        printf("P%d i=%d\n", proc_num, i);
    }*/



