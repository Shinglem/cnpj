package io.github.shinglem.cnpj

import io.github.shinglem.cnpj.webserver.WebserverVerticle
import io.vertx.core.DeploymentOptions
import io.vertx.core.Vertx
import org.slf4j.LoggerFactory


class App {

    companion object {
        final val logger = LoggerFactory.getLogger(this::class.java.name)
    }


}


fun main(args: Array<String>) {
    System.setProperty("vertx.logger-delegate-factory-class-name", "io.vertx.core.logging.SLF4JLogDelegateFactory")
    System.setProperty("user.timezone", "GMT +08")
    System.setProperty("kotlinx.coroutines.debug", "")
    System.setProperty("vertx.disableFileCaching", "false")

    val logger = App.logger

    val vertx = VertxInstance.vertx


    vertx.deployVerticle(WebserverVerticle::class.java , DeploymentOptions()){
        if (it.succeeded()) {
            logger.info("server start")
        }else{
            logger.error("error !!!")
        }
    }



}

object VertxInstance {
    val vertx = Vertx.vertx()
}