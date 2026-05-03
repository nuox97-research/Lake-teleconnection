function Q = modularity(A, communities)
    % A: 加权邻接矩阵
    % communities: 节点的社区分配（一个向量，表示每个节点所属的社区）

    % 计算总边权重 m
    m = sum(A(:)) / 2;  % 所有边的权重之和

    % 初始化模块度 Q
    Q = 0;

    % 获取节点数量
    N = size(A, 1);

    % 计算模块度
    for i = 1:N
        for j = 1:N
            if communities(i) == communities(j)
                % 节点 i 和节点 j 在同一社区时，累加模块度
                Q = Q + A(i, j) - (sum(A(i, :)) * sum(A(j, :))) / (2 * m);
            end
        end
    end

    % 最终的模块度值
    Q = Q / (2 * m);
end